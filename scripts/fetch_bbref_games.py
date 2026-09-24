#!/usr/bin/env python3
"""
Download NBA regular-season results from basketball-reference.com and write
them to data/gamesYYZZ.csv in the format used by the Tanking package, i.e.,

    Date,Time,Team1,Score1,Team2,Score2,Label,OT,Attendance,Misc

where Team1 is the visiting team and Team2 is the home team.

Only regular-season games are kept:
  * everything after the "Playoffs" separator row is dropped;
  * play-in games are dropped;
  * the NBA Cup (In-Season Tournament) championship game, which does not
    count in the regular-season standings (2023-24 onwards), is dropped.
Each season is validated (every team plays 82 games, 1230 games total)
before it is written.

After downloading, data/winpct.csv (win pct of the team in each final
position, one column per season) is regenerated from all data/games*.csv.

Usage (from the project root; needs only the Python standard library):
    python3 scripts/fetch_bbref_games.py 2022 2023 2024 2025 2026
    python3 scripts/fetch_bbref_games.py --winpct-only

A season is named by the year in which it ends (2022 = 2021-22).
basketball-reference rate-limits scrapers (~20 requests/minute), so requests
are spaced out by --delay seconds.
"""

import argparse
import collections
import csv
import datetime
import glob
import html.parser
import os
import re
import sys
import time
import urllib.request

BASE_URL = "https://www.basketball-reference.com"
HEADER = ["Date", "Time", "Team1", "Score1", "Team2", "Score2",
          "Label", "OT", "Attendance", "Misc"]
NUM_TEAMS = 30
NUM_TEAM_GAMES = 82

# NBA Cup championship games (not counted in regular-season standings);
# used only if the game cannot be identified from the Notes column
CUP_FINAL_DATES = {
    2024: datetime.date(2023, 12, 9),
    2025: datetime.date(2024, 12, 17),
    2026: datetime.date(2025, 12, 16),
}

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_DIR, "data")


def fetch(url, delay):
    time.sleep(delay)
    print(f"Fetching {url}", file=sys.stderr)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (tanking research script)"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


class ScheduleParser(html.parser.HTMLParser):
    """Collect rows of <table id="schedule"> as dicts keyed by data-stat."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.in_table = False
        self.table_depth = 0
        self.row = None
        self.row_is_thead = False
        self.cell_stat = None
        self.cell_text = []
        self.rows = []  # list of dicts, or the string "SEPARATOR:<text>"

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "table":
            if self.in_table:
                self.table_depth += 1
            elif attrs.get("id") == "schedule":
                self.in_table = True
                self.table_depth = 1
            return
        if not self.in_table:
            return
        if tag == "tr":
            self.row = {}
            self.row_is_thead = "thead" in (attrs.get("class") or "").split()
        elif tag in ("td", "th") and self.row is not None:
            self.cell_stat = attrs.get("data-stat") or ""
            self.cell_text = []

    def handle_endtag(self, tag):
        if not self.in_table:
            return
        if tag == "table":
            self.table_depth -= 1
            if self.table_depth == 0:
                self.in_table = False
        elif tag in ("td", "th") and self.row is not None and self.cell_stat is not None:
            self.row[self.cell_stat] = " ".join("".join(self.cell_text).split())
            self.cell_stat = None
        elif tag == "tr" and self.row is not None:
            if self.row_is_thead:
                self.rows.append("SEPARATOR:" + " ".join(self.row.values()))
            elif self.row:
                self.rows.append(self.row)
            self.row = None

    def handle_data(self, data):
        if self.in_table and self.cell_stat is not None:
            self.cell_text.append(data)


def get(row, *keys):
    for k in keys:
        if k in row:
            return row[k]
    raise KeyError(f"None of {keys} found in row {row}")


def parse_date(s):
    # e.g. "Tue, Oct 19, 2021"
    return datetime.datetime.strptime(s.replace(",", ""), "%a %b %d %Y").date()


def parse_schedule_page(page):
    # Some bbref tables are wrapped in HTML comments; unwrap them just in case
    page = page.replace("<!--", "").replace("-->", "")
    p = ScheduleParser()
    p.feed(page)
    games = []
    saw_playoffs = False
    for row in p.rows:
        if isinstance(row, str):
            if "playoffs" in row.lower():
                saw_playoffs = True
                break  # everything after this is postseason
            continue
        if "date_game" not in row and "date" not in row:
            continue  # header row repeated inside tbody
        date_str = get(row, "date_game", "date")
        if not date_str or date_str.lower() == "date":
            continue
        if date_str.lower() == "playoffs":
            saw_playoffs = True
            break
        games.append({
            "date": parse_date(date_str),
            "date_str": date_str.replace(",", ""),
            "time": row.get("game_start_time", ""),
            "visitor": get(row, "visitor_team_name"),
            "visitor_pts": get(row, "visitor_pts"),
            "home": get(row, "home_team_name"),
            "home_pts": get(row, "home_pts"),
            "label": row.get("box_score_text", ""),
            "ot": row.get("overtimes", ""),
            "attendance": row.get("attendance", "").replace(",", ""),
            "notes": row.get("game_remarks", ""),
        })
    return games, saw_playoffs


def month_urls(season, delay):
    index_url = f"{BASE_URL}/leagues/NBA_{season}_games.html"
    page = fetch(index_url, delay)
    links = re.findall(rf'href="(/leagues/NBA_{season}_games-[a-z]+\.html)"', page)
    urls = []
    for link in links:
        if link not in urls:
            urls.append(link)
    if not urls:
        # Index page itself holds the first month; no month filter found
        return [index_url], page
    return [BASE_URL + u for u in urls], None


def fetch_season(season, delay):
    urls, first_page = month_urls(season, delay)
    games = []
    for url in urls:
        page = first_page if (first_page is not None and url.endswith("_games.html")) else fetch(url, delay)
        month_games, saw_playoffs = parse_schedule_page(page)
        games.extend(month_games)
        # The playoffs start in the month in which the regular season ends;
        # later months only have playoff games, so stop once we see them.
        if saw_playoffs:
            break
    return games


def filter_regular_season(season, games):
    """Remove play-in games and the NBA Cup final; validate the result."""
    kept = []
    for g in games:
        if not g["visitor_pts"] or not g["home_pts"]:
            print(f"  [{season}] skipping unplayed game: {g['date_str']} {g['visitor']} @ {g['home']}", file=sys.stderr)
            continue
        if re.search(r"play-?in", g["notes"], re.I):
            print(f"  [{season}] dropping play-in game: {g['date_str']} {g['visitor']} @ {g['home']} ({g['notes']})", file=sys.stderr)
            continue
        kept.append(g)

    # NBA Cup championship game
    cup_final = [g for g in kept
                 if re.search(r"(cup|tournament)", g["notes"], re.I)
                 and re.search(r"(final|championship)", g["notes"], re.I)
                 and not re.search(r"(semi|quarter)", g["notes"], re.I)]
    if not cup_final and season in CUP_FINAL_DATES:
        cup_final = [g for g in kept if g["date"] == CUP_FINAL_DATES[season]]
        if len(cup_final) != 1:
            sys.exit(f"[{season}] expected exactly one game on NBA Cup final date {CUP_FINAL_DATES[season]}, found {len(cup_final)}")
    for g in cup_final:
        print(f"  [{season}] dropping NBA Cup final: {g['date_str']} {g['visitor']} @ {g['home']} ({g['notes']})", file=sys.stderr)
    kept = [g for g in kept if g not in cup_final]

    # Anything a team plays after its 82nd game (e.g. unlabeled play-in games) is not regular season
    kept.sort(key=lambda g: g["date"])  # stable: keeps bbref order within a day
    count = collections.Counter()
    last_regular_date = None
    regular, extra = [], []
    for g in kept:
        if count[g["visitor"]] >= NUM_TEAM_GAMES or count[g["home"]] >= NUM_TEAM_GAMES:
            extra.append(g)
            continue
        count[g["visitor"]] += 1
        count[g["home"]] += 1
        last_regular_date = g["date"]
        regular.append(g)
    for g in extra:
        if g["date"] < last_regular_date:
            sys.exit(f"[{season}] game after a team's 82nd game occurs before the end of the regular season: {g}")
        print(f"  [{season}] dropping post-regular-season game: {g['date_str']} {g['visitor']} @ {g['home']} ({g['notes']})", file=sys.stderr)

    teams = set(count)
    if len(teams) != NUM_TEAMS:
        sys.exit(f"[{season}] found {len(teams)} teams, expected {NUM_TEAMS}: {sorted(teams)}")
    bad = {t: c for t, c in count.items() if c != NUM_TEAM_GAMES}
    if bad:
        sys.exit(f"[{season}] teams without {NUM_TEAM_GAMES} games: {bad}")
    assert len(regular) == NUM_TEAMS * NUM_TEAM_GAMES // 2
    for g in regular:
        if int(g["visitor_pts"]) == int(g["home_pts"]):
            sys.exit(f"[{season}] tied game: {g}")
    return regular


def season_filename(season):
    return os.path.join(DATA_DIR, f"games{(season - 1) % 100:02d}{season % 100:02d}.csv")


def write_season(season, games):
    fname = season_filename(season)
    with open(fname, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        for g in games:
            w.writerow([g["date_str"], g["time"], g["visitor"], g["visitor_pts"],
                        g["home"], g["home_pts"], g["label"], g["ot"],
                        g["attendance"], g["notes"].replace(",", ";")])
    print(f"Wrote {len(games)} games to {fname}", file=sys.stderr)


def season_start_year(fname):
    yy = int(os.path.basename(fname)[5:7])
    return yy + (1900 if yy > 90 else 2000)


def write_winpct():
    """Regenerate data/winpct.csv: column = season (newest first), row = final position."""
    files = sorted(glob.glob(os.path.join(DATA_DIR, "games[0-9][0-9][0-9][0-9].csv")),
                   key=season_start_year, reverse=True)
    header, columns = [], []
    for fname in files:
        header.append(f"{season_start_year(fname)}-{os.path.basename(fname)[7:9]}")
        wins, played = collections.Counter(), collections.Counter()
        with open(fname, newline="") as f:
            for row in list(csv.reader(f))[1:]:
                s1, s2 = int(row[3]), int(row[5])
                if s1 == s2:  # canceled game (e.g., 2012-13 BOS-IND)
                    continue
                played[row[2]] += 1
                played[row[4]] += 1
                wins[row[2] if s1 > s2 else row[4]] += 1
        columns.append(sorted((100.0 * wins[t] / played[t] for t in played), reverse=True))
    fname = os.path.join(DATA_DIR, "winpct.csv")
    with open(fname, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        for pos in range(NUM_TEAMS):
            w.writerow([f"{col[pos]:.8f}" for col in columns])
    print(f"Wrote {fname} with seasons {header}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seasons", nargs="*", type=int, help="season end years, e.g. 2022 for 2021-22")
    ap.add_argument("--delay", type=float, default=4.0, help="seconds between requests (default 4)")
    ap.add_argument("--winpct-only", action="store_true", help="only regenerate data/winpct.csv")
    args = ap.parse_args()

    if not args.winpct_only:
        if not args.seasons:
            ap.error("give at least one season, or --winpct-only")
        for season in args.seasons:
            games = fetch_season(season, args.delay)
            write_season(season, filter_regular_season(season, games))
    write_winpct()


if __name__ == "__main__":
    main()
