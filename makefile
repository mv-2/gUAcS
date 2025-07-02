# Prepares workspace for running
build_install:
	pip install .

# maturin build after code updates
build:
	maturin develop --release

# run project and all required build steps
run:
	make build
	python examples/func_driver.py

# clean all generated files and build
clean_build:
	cargo clean
	make build

# run from clean
clean_run:
	make clean_build
	python src/driver.py

# build and install all packages required for working remotely on new machines
install_all:
	# install rust and required components for development
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
	rustup component add rust-analyzer
	rustup component add clippy
	rustup component add rustfmt
	rustup component add cargo
	
	# install python stuff and build project
	make build_install
