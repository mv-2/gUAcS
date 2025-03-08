# writes config file
cfg_write:
	python src/make_config.py

# Prepares workspace for running
build:
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
