include ../default_targets.h

# run this test only if a grid producing executable exists
test:
	@(if [ -x "$(DELPHI_PATH)" ] ; then \
	        echo "This DelPhi executable was found:" ;\
	        echo "    $(DELPHI_PATH)" ;\
	        make chemgrid_generation ;\
	        make delphi_grid_generation ;\
	        cp -p gconv_delphi.in gconv.in;\
	        make grid_convert ;\
	elif [ -x "$(QNIFFT_PATH)" ] ; then \
	        echo "This Qnifft executable was found:" ;\
	        echo "    $(QNIFFT_PATH)" ;\
	        make chemgrid_generation ;\
	        make qnifft_grid_generation ;\
	        cp -p gconv_qnifft.in gconv.in;\
	        make grid_convert ;\
	else \
	        echo "" ;\
	        echo "An executable to produce grids was not found." ;\
	        echo "    DELPHI_PATH = $(DELPHI_PATH)" ;\
	        echo "    QNIFFT_PATH = $(QNIFFT_PATH)" ;\
	        echo "    Define DELPHI_PATH or QNIFFT_PATH to its full path." ;\
	        echo "    Skipping this test." ;\
	        echo "" ;\
	fi ;\
	)


clean: clean_chemgrid_generation clean_delphi_grid_generation clean_grid_convert
	/bin/rm -f gconv.in

