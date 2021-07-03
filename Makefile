JULIA := julia --color=yes --project=@.

build/sysimage.so: build/sysimage.jl Project.toml Manifest.toml
	mkdir -p build; mkdir -p results/tmp
	$(JULIA) --trace-compile="build/precompile.jl" test/test_script.jl 1
	$(JULIA) build/sysimage.jl
