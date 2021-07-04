JULIA := julia --color=yes --project=@.
SYSIMAGE = build/JuliaTanking.so

$(SYSIMAGE): build/sysimage.jl Project.toml Manifest.toml
	mkdir -p build; mkdir -p results/tmp
	$(JULIA) --trace-compile="build/precompile.jl" test/test_script.jl 1
	$(JULIA) build/sysimage.jl

test: $(SYSIMAGE)
	@echo Testing Tanking package through sysimage
	$(JULIA) --sysimage $(SYSIMAGE) test/test_script.jl 1

.PHONY: test
