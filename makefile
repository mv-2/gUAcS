# writes config file
cfg_write:
	python src/make_config.py

# Prepares workspace for running
build:
	make install_all
	make cfg_write
	maturin develop --release

# run project and all required build steps
run:
	make build
	python src/driver.py

# clean all generated files and build
clean_build:
	cargo clean
	make build

# run from clean
clean_run:
	make clean_build
	python src/driver.py

# build and install all dependancies on new system
# TODO: Make generalised build and install script for all system types
install_all:
	# install rust and required components for development
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
	rustup component add rust-analyzer
	rustup component add clippy
	rustup component add rustfmt
	rustup component add cargo
	
	# install all python and required components
	mkdir -p ~/miniconda3
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
	bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
	rm ~/miniconda3/miniconda.sh
	source ~/miniconda3/bin/activate
	conda init --all
	conda create --name guac_env --file py_req.txt
	conda activate guac_env
	
