SHELL=/bin/bash
DIR := ${CURDIR}
LIB = $${LD_LIBRARY_PATH}
RDKITLIB = $${LD_LIBRARY_PATH}

install: rdkit

getdepend:
	wget -c -nv https://archives.boost.io/release/1.71.0/source/boost_1_71_0.tar.gz
	wget -c -nv https://github.com/rdkit/rdkit/archive/Release_2019_09_1.tar.gz
	wget -c -nv https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
	mv Release_2019_09_1.tar.gz rdkit-Release_2019_09_1.tar.gz

conda:
	conda create -n dock6_py37 python=3.7 -yq
	source ${CONDABASE}/bin/activate dock6_py37 && \
	    ln -s ${CONDABASE}/envs/dock6_py37/include/python3.7m \
	          ${CONDABASE}/envs/dock6_py37/include/python3.7  && \
	    ${CONDABASE}/bin/conda install -yq numpy && \
	    ${CONDABASE}/bin/conda install -yq cmake=3.11.1
	
extract:
	tar -xf rdkit-Release_2019_09_1.tar.gz 
	tar -xf boost_1_71_0.tar.gz
	tar -xf eigen-3.3.9.tar.gz

boostcompile:
	source ${CONDABASE}/bin/activate dock6_py37 && cd ./boost_1_71_0/ && \
	    ./bootstrap.sh --prefix=../boost
	source ${CONDABASE}/bin/activate dock6_py37 && cd ./boost_1_71_0/ && \
	    ./b2 install --prefix=../boost

rdkitcompile:
	source ${CONDABASE}/bin/activate dock6_py37 && \
	    cd rdkit-Release_2019_09_1 && \
	    rm -rf build && mkdir build && cd build && \
	    cmake .. -DRDK_INSTALL_INTREE=OFF \
	             -DCMAKE_INSTALL_PREFIX=${DIR}/rdkit-Release_2019_09_1 \
	             -DBOOST_ROOT=${DIR}/boost \
	             -DBoost_NO_SYSTEM_PATHS=ON \
	             -DBoost_USE_STATIC_LIBS=OFF \
	             -DEIGEN3_INCLUDE_DIR=${DIR}/eigen-3.3.9 \
	             -DPYTHON_EXECUTABLE=${CONDABASE}/envs/dock6_py37/bin/python3 \
	             -DPYTHON_NUMPY_INCLUDE_PATH=${CONDABASE}/envs/dock6_py37/lib/python3.7/site-packages/numpy/core/include \
	             -DRDK_BUILD_PYTHON_WRAPPERS=OFF && \
	    make && \
	    make install


echobegin:
	@echo " "
	@echo "RDKIT COMPILATION STARTING..."
	@echo " "

echoend:
	@echo " "
	@echo "######################## RDKIT COMPILATION COMPLETE" 
	@echo "######################## on `date` "
	@echo "######################## Please export these variables below in your bashrc file"
	@echo " "
	@echo "    export RDBASE=${DIR}/rdkit-Release_2019_09_1"
	@echo '    export LD_LIBRARY_PATH=${LIB}:${DIR}/rdkit-Release_2019_09_1/lib' 
	@echo '    export LD_LIBRARY_PATH=${DIR}/boost/lib:${LIB}'
	@echo "    export BOOST=${DIR}/boost"
	@echo " "
	@echo "############ You must export these variables or DOCK will not compile with RDKIT ##########"
	@echo " "


rdkitclean:
	/bin/rm -f boost_1_71_0.tar.gz eigen-3.3.9.tar.gz Release_2019_09_1.tar.gz \
	           rdkit-Release_2019_09_1.tar.gz
	/bin/rm -rf ./boost ./boost_1_71_0 ./eigen-3.3.9 ./rdkit-Release_2019_09_1

rdkit:
	@echo "RDKIT COMPILATION CAN TAKE HOURS"
	@echo "IF YOU HAVE YOUR OWN COMPILED VERSION PLEASE EXPORT YOUR VARIABLES IN YOUR BASHRC OR EQUIVALENT FILE"
	@echo "EXPORT VARIABLE TEMPLATE SHOWN BELOW"
	@echo ""
	@echo "    export RDBASE=PATH/TO/RDKIT/ROOT/DIR"
	@echo '    export LD_LIBRARY_PATH=${RDKITLIB}:PATH/TO/RDKIT/ROOT/DIR/lib' 
	@echo '    export LD_LIBRARY_PATH=PATH/TO/BOOST/ROOT/DIR/lib:${RDKITLIB}'
	@echo "    export BOOST=PATH/TO/BOOST/ROOT/DIR/"
	@echo ""
	@echo "############################################################################"
	@echo "BEFORE RUNNING RDKIT COMPILATION, YOU MUST BE IN A CONDA ENVIRONMENT."
	@echo "FURTHER,YOU MUST EXPORT ANACONDA TO CONDABASE IN YOUR BASHRC OR EQUIVALENT FILE. SHOWN BELOW."
	@echo "" 
	@echo "    export CONDABASE=PATH/TO/CONDA/ROOT/"
	@echo " " 
	@echo "RUN RDKIT COMPILATION [Y/n]?" 
	@read rdkitconfirm; \
	if [ "$${rdkitconfirm}" = "Y" ]; then \
	    $(MAKE) echobegin; \
	    $(MAKE) getdepend; \
	    $(MAKE) conda; \
	    $(MAKE) extract; \
	    $(MAKE) boostcompile; \
	    $(MAKE) rdkitcompile; \
	    $(MAKE) echoend; \
	else \
	    echo " you chose to exit. Exiting..."; \
	fi

