ifneq ($(findstring x86_64-linux,$(tgt_arch)),)
qtincdir  := qt/include_64
else
qtincdir  := qt/include
endif

# List targets (if any) for this package
tgtnames := online_ami offline_ami blviewer bldipimbclient qttest

# List source files for each target
tgtsrcs_online_ami += qtclient.cc

tgtsrcs_offline_ami := offline_ami.cc XtcFileClient.cc XtcFileClient_moc.cc

tgtsrcs_blviewer := blvclient.cc blvclient_moc.cc

tgtsrcs_bldipimbclient := bldIpimbClient.cc bldIpimbClient_moc.cc

tgtsrcs_qttest := qttest.cc

# List system libraries (if any) needed by exe_a as <dir>/<lib>. 
# Note that <lib> is the name of the library, not of the file: i.e.
# <lib> for 'libc.so' is 'c'. Low level first.
# tgtslib_exe_a := $(USRLIBDIR)/rt

qt_libs := qt/QtGui qt/QtCore
qt_libs += qwt/qwt

# List project libraries (if any) needed by exe_a as <project>/<lib>.
# Note that <lib> is the name of the library, not of the file: i.e.
# <lib> for 'libc.so' is 'c'. Low level first.
tgtlibs_online_ami := pdsdata/xtcdata pdsdata/acqdata
tgtlibs_online_ami += ami/service ami/data ami/server ami/client ami/calib ami/amiqt
tgtlibs_online_ami += $(qt_libs)

datalibs := pdsdata/xtcdata pdsdata/opal1kdata pdsdata/quartzdata pdsdata/pulnixdata pdsdata/camdata pdsdata/pnccddata pdsdata/evrdata pdsdata/acqdata pdsdata/controldata pdsdata/princetondata pdsdata/ipimbdata pdsdata/encoderdata pdsdata/fccddata pdsdata/lusidata pdsdata/cspaddata pdsdata/xampsdata pdsdata/fexampdata pdsdata/gsc16aidata pdsdata/epics pdsdata/usdusbdata
datalibs += pdsdata/timepixdata
datalibs += pdsdata/phasicsdata
datalibs += pdsdata/cspad2x2data
datalibs += pdsdata/oceanopticsdata pdsdata/flidata pdsdata/andordata pdsdata/orcadata

#
# Need all pdsdata libraries to support dynamic linking of plug-in modules
#
tgtlibs_offline_ami := $(datalibs) pdsdata/appdata pdsdata/anadata pdsdata/indexdata
tgtlibs_offline_ami += ami/service ami/data ami/server ami/client ami/calib ami/event ami/app ami/amiqt
tgtlibs_offline_ami += $(qt_libs)

tgtlibs_blviewer := $(datalibs) pdsapp/configdb pdsapp/configdbg
tgtlibs_blviewer += ami/service ami/data ami/server ami/client ami/calib ami/amiqt
tgtlibs_blviewer += $(qt_libs)

tgtlibs_bldipimbclient := $(datalibs) pdsapp/configdb pdsapp/configdbg
tgtlibs_bldipimbclient += ami/service ami/data ami/server ami/client ami/calib ami/amiqt
tgtlibs_bldipimbclient += $(qt_libs)

tgtlibs_qttest := $(datalibs)
tgtlibs_qttest += ami/service ami/data ami/server ami/client ami/calib ami/amiqt
tgtlibs_qttest += $(qt_libs)

# List special include directories (if any) needed by exe_a as
# <project>/<incdir>. Note that the top level release directory is
# already in the search path.
qt_incs := $(qtincdir) qwt/include
# qwt includes qt headers without package prefix!
qt_incs += $(qtincdir)/Qt
qt_incs += $(qtincdir)/Qt3Support
qt_incs += $(qtincdir)/QtAssistant
qt_incs += $(qtincdir)/QtCore
qt_incs += $(qtincdir)/QtDesigner
qt_incs += $(qtincdir)/QtGui
qt_incs += $(qtincdir)/QtNetwork
qt_incs += $(qtincdir)/QtOpenGL
qt_incs += $(qtincdir)/QtScript
qt_incs += $(qtincdir)/QtSql
qt_incs += $(qtincdir)/QtSvg
qt_incs += $(qtincdir)/QtTest
qt_incs += $(qtincdir)/QtUiTools
qt_incs += $(qtincdir)/QtXml

tgtincs_online_ami  := $(qt_incs)
tgtincs_offline_ami := $(qt_incs)
tgtincs_blviewer    := $(qt_incs)
tgtincs_bldipimbclient   := $(qt_incs)
tgtincs_qttest := $(qt_incs)

# List system include directories (if any) needed by exe_a as <incdir>.
# tgtsinc_exe_a := /usr/include

# List libraries (if any) for this package
libnames := amiqt

# List source files for each library
libsrcs_amiqt := Status.cc Status_moc.cc
libsrcs_amiqt += Path.cc
libsrcs_amiqt += FeatureRegistry.cc FeatureRegistry_moc.cc
libsrcs_amiqt += SMPRegistry.cc SMPRegistry_moc.cc
libsrcs_amiqt += QHComboBox.cc QHComboBox_moc.cc
libsrcs_amiqt += FeatureBox.cc FeatureBox_moc.cc
libsrcs_amiqt += QtTree.cc QtTree_moc.cc
libsrcs_amiqt += RunTree.cc RunTree_moc.cc
libsrcs_amiqt += FeatureTree.cc FeatureTree_moc.cc
libsrcs_amiqt += CExpression.cc
libsrcs_amiqt += Calculator.cc Calculator_moc.cc
libsrcs_amiqt += FeatureCalculator.cc
libsrcs_amiqt += DescTH2T.cc DescTH2T_moc.cc
libsrcs_amiqt += DescTH2F.cc DescTH2F_moc.cc
libsrcs_amiqt += DescTH1F.cc DescTH1F_moc.cc
libsrcs_amiqt += DescBinning.cc DescBinning_moc.cc
libsrcs_amiqt += DescProf.cc DescProf_moc.cc
libsrcs_amiqt += DescScan.cc DescScan_moc.cc
libsrcs_amiqt += DescChart.cc DescChart_moc.cc
libsrcs_amiqt += ScalarPlotDesc.cc
libsrcs_amiqt += XYHistogramPlotDesc.cc
libsrcs_amiqt += XYProjectionPlotDesc.cc
libsrcs_amiqt += RPhiProjectionPlotDesc.cc
libsrcs_amiqt += QtPersistent.cc QtUtils.cc
libsrcs_amiqt += ChannelDefinition.cc ChannelDefinition_moc.cc
libsrcs_amiqt += AxisArray.cc
libsrcs_amiqt += AxisBins.cc
libsrcs_amiqt += AxisControl.cc AxisControl_moc.cc
libsrcs_amiqt += Control.cc Control_moc.cc
libsrcs_amiqt += TransformConstant.cc TransformConstant_moc.cc
libsrcs_amiqt += Transform.cc Transform_moc.cc
libsrcs_amiqt += Condition.cc Condition_moc.cc
libsrcs_amiqt += Contour.cc
libsrcs_amiqt += Filter.cc Filter_moc.cc
libsrcs_amiqt += ChannelMath.cc ChannelMath_moc.cc
libsrcs_amiqt += QtBase.cc
libsrcs_amiqt += QtPlotCurve.cc QtPlotCurve_moc.cc
libsrcs_amiqt += QtWaveform.cc QtTH2F.cc QtTH1F.cc QtProf.cc QtChart.cc QtImage.cc QtScan.cc QtChannelMath.cc QtEmpty.cc QtTable.cc QtTable_moc.cc
libsrcs_amiqt += PlotFactory.cc
libsrcs_amiqt += PlotFrame.cc PlotFrame_moc.cc
libsrcs_amiqt += PrintAction.cc PrintAction_moc.cc
libsrcs_amiqt += WaveformDisplay.cc WaveformDisplay_moc.cc
libsrcs_amiqt += CursorDefinition.cc CursorDefinition_moc.cc
libsrcs_amiqt += QtPlot.cc QtPlot_moc.cc
libsrcs_amiqt += QtPlotStyle.cc
libsrcs_amiqt += QtPlotSelector.cc QtPlotSelector_moc.cc
libsrcs_amiqt += QtOverlay.cc QtOverlay_moc.cc
libsrcs_amiqt += CursorPlot.cc CursorPlot_moc.cc
libsrcs_amiqt += CursorPost.cc CursorOverlay.cc
libsrcs_amiqt += CursorsX.cc CursorsX_moc.cc
libsrcs_amiqt += EdgePlot.cc EdgePlot_moc.cc
libsrcs_amiqt += EdgeOverlay.cc
libsrcs_amiqt += EdgeCursor.cc EdgeCursor_moc.cc
libsrcs_amiqt += EdgeFinder.cc EdgeFinder_moc.cc
libsrcs_amiqt += CurveFitPost.cc CurveFitOverlay.cc
libsrcs_amiqt += CurveFitPlot.cc CurveFitPlot_moc.cc
libsrcs_amiqt += CurveFit.cc CurveFit_moc.cc
libsrcs_amiqt += QtPWidget.cc QtPWidget_moc.cc
libsrcs_amiqt += AbsClient.cc AbsClient_moc.cc
libsrcs_amiqt += Client.cc Client_moc.cc
libsrcs_amiqt += WaveformClient.cc
libsrcs_amiqt += ZoomPlot.cc ZoomPlot_moc.cc
libsrcs_amiqt += PeakFit.cc PeakFit_moc.cc
libsrcs_amiqt += PeakFitPlot.cc PeakFitPlot_moc.cc
libsrcs_amiqt += PeakFitPost.cc PeakFitOverlay.cc
libsrcs_amiqt += PeakPlot.cc PeakPlot_moc.cc
libsrcs_amiqt += PeakFinder.cc PeakFinder_moc.cc
libsrcs_amiqt += BlobFinder.cc BlobFinder_moc.cc
libsrcs_amiqt += ProjectionPlot.cc ProjectionPlot_moc.cc
libsrcs_amiqt += TwoDPlot.cc TwoDPlot_moc.cc
libsrcs_amiqt += ImageGrid.cc
libsrcs_amiqt += ImageFrame.cc ImageFrame_moc.cc
libsrcs_amiqt += ImageColorControl.cc ImageColorControl_moc.cc
libsrcs_amiqt += ImageXYControl.cc ImageXYControl_moc.cc
libsrcs_amiqt += ColorMaps.cc
libsrcs_amiqt += ImageScale.cc ImageScale_moc.cc
libsrcs_amiqt += ImageDisplay.cc ImageDisplay_moc.cc
libsrcs_amiqt += Rect.cc RectROI.cc RectROI_moc.cc
libsrcs_amiqt += RectangleCursors.cc RectangleCursors_moc.cc
libsrcs_amiqt += AnnulusCursors.cc AnnulusCursors_moc.cc
libsrcs_amiqt += ImageIntegral.cc ImageContrast.cc
libsrcs_amiqt += ImageFunctions.cc ImageFunctions_moc.cc
libsrcs_amiqt += ImageContourProjection.cc ImageContourProjection_moc.cc
libsrcs_amiqt += ImageRPhiProjection.cc ImageRPhiProjection_moc.cc
libsrcs_amiqt += ImageXYProjection.cc ImageXYProjection_moc.cc
libsrcs_amiqt += CrossHair.cc CrossHair_moc.cc
libsrcs_amiqt += CrossHairDelta.cc CrossHairDelta_moc.cc
libsrcs_amiqt += ImageGridScale.cc ImageGridScale_moc.cc
libsrcs_amiqt += ImageClient.cc
libsrcs_amiqt += CspadClient.cc CspadClient_moc.cc
libsrcs_amiqt += FccdClient.cc FccdClient_moc.cc
libsrcs_amiqt += EnvPlot.cc EnvPlot_moc.cc EnvTable.cc EnvTable_moc.cc
libsrcs_amiqt += EnvPost.cc EnvOverlay.cc
libsrcs_amiqt += PostAnalysis.cc
libsrcs_amiqt += EnvClient.cc EnvClient_moc.cc
libsrcs_amiqt += TdcPlot.cc TdcPlot_moc.cc
libsrcs_amiqt += TdcClient.cc TdcClient_moc.cc
libsrcs_amiqt += SummaryClient.cc SummaryClient_moc.cc
#libsrcs_amiqt += DetectorGroup.cc DetectorGroup_moc.cc
#libsrcs_amiqt += DetectorSave.cc
#libsrcs_amiqt += DetectorReset.cc
libsrcs_amiqt += DetectorListItem.cc
libsrcs_amiqt += Defaults.cc
libsrcs_amiqt += RateCalculator.cc RateCalculator_moc.cc
libsrcs_amiqt += RateDisplay.cc RunMaster.cc
libsrcs_amiqt += FilterSetup.cc FilterSetup_moc.cc
libsrcs_amiqt += PWidgetManager.cc PWidgetManager_moc.cc
libsrcs_amiqt += DetectorSelect.cc DetectorSelect_moc.cc
libsrcs_amiqt += MaskFrame.cc MaskFrame_moc.cc
libsrcs_amiqt += MaskDisplay.cc MaskDisplay_moc.cc

libsrcs_amiblv := blv
# List special include directories (if any) needed by lib_a as
# <project>/<incdir>. Note that the top level release directory is
# already in the search path.
libincs_amiqt := $(qt_incs)

# List system include directories (if any) needed by lib_a as <incdir>.
# libsinc_lib_a := /usr/include


