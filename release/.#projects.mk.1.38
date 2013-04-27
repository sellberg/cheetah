# List of projects (low level first)

#
# 32-bit linux
#
ifneq ($(findstring i386-linux,$(tgt_arch)),)
projects := pdsdata \
      acqiris \
      evgr \
      leutron \
      edt \
      qt \
      qwt \
      epics \
      offlinedb \
      pvcam \
      relaxd \
      fli \
      andor \
      libusb \
      usdusb4

epics_use   := /reg/g/pcds/package/external/epicsca-pcds-R1.0-r410
acqiris_use := /reg/g/pcds/package/external/acqiris_3.3a
evgr_use := /reg/g/pcds/package/external/evgr_V00-00-02
leutron_use := /reg/g/pcds/package/external/leutron_V00-00-00
edt_use := /reg/g/pcds/package/external/edt
relaxd_use := /reg/g/pcds/package/external/relaxd-1.8.0
qt_use := /reg/g/pcds/package/external/qt-4.3.4
qwt_use := /reg/g/pcds/package/external/qwt-5.1.1-wfopt
offlinedb_use := /reg/g/pcds/package/external/offlinedb-1.5.1
pvcam_use := /reg/g/pcds/package/external/pvcam2.7.1.7
fli_use   := /reg/g/pcds/package/external/fli-dist-1.71
andor_use := /reg/g/pcds/package/external/andor-2.93.30007
libusb_use := /reg/g/pcds/package/external/libusb-1.0.0
usdusb4_use := /reg/g/pcds/package/external/usdusb4

pdsdata_use := release

endif

#
# 64-bit linux
#
ifneq ($(findstring x86_64-linux,$(tgt_arch)),)
projects := pdsdata \
      qt \
      qwt \
      epics \
      evgr \
      edt \
      offlinedb \
      leutron \
      python \
      libraw1394 \
      libdc1394 \
      fli \
      andor \
      libusb \
      usdusb4

epics_use   := /reg/g/pcds/package/external/epicsca-pcds-R1.0-r410
evgr_use    := /reg/g/pcds/package/external/evgr_V00-00-02
qt_use      := /reg/g/pcds/package/external/qt-4.3.4
qwt_use     := /reg/g/pcds/package/external/qwt-5.1.1-wfopt
python_use  := /reg/g/pcds/package/python-2.5.2
libraw1394_use := /reg/g/pcds/package/external/libdc1394
libdc1394_use := /reg/g/pcds/package/external/libdc1394
offlinedb_use := /reg/g/pcds/package/external/offlinedb-1.5.1
edt_use := /reg/g/pcds/package/external/edt
leutron_use := /reg/g/pcds/package/external/leutron_V00-00-00
fli_use   := /reg/g/pcds/package/external/fli-dist-1.71
andor_use := /reg/g/pcds/package/external/andor-2.93.30007
libusb_use := /reg/g/pcds/package/external/libusb-1.0.0
usdusb4_use := /reg/g/pcds/package/external/usdusb4

pds_use := release
pdsdata_use := release

endif
