# List of packages (low level first)
ifneq ($(findstring ppc-rtems-rce,$(tgt_arch)),)
packages := xtc
else
packages := ipimb encoder pnCCD acqiris camera evr opal1k pulnix control xtc \
            epics bld princeton fccd cspad lusi xamps fexamp gsc16ai \
            phasics timepix cspad2x2 oceanoptics fli quartz andor usdusb orca \
	    compress ana index app
endif
