#!MC 1410
$!OpenLayout  "/home/geve/Dropbox/EnKF_advection/Run/review.lay"

$!VARSET |VAR|=0
$!LOOP 10
$!VARSET |VAR| += 10
$!ReadDataSet  '"/home/geve/Dropbox/EnKF_advection/Run/Solution/sol|VAR|.00A.dat" "/home/geve/Dropbox/EnKF_advection/Run/Solution/obs|VAR|.00A.dat" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"V1" "V2" "V3" "V4" "V5"'
$!AlterData 
  Equation = 'V6=2.0*V5'
$!LinePlotLayers ShowErrorBars = Yes
$!LineMap [1]  ErrorBars{Var = 6}
$!LineMap [5]  ErrorBars{Var = 6}
$!AttachText 
  AnchorPos
    {
    X = 9.44
    Y = 89.65
    }
  TextShape
    {
    IsBold = No
    SizeUnits = Frame
    Height = 6
    }
  Text = 'T=|VAR|'
$!Redraw 
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup ImageWidth = 1569
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ExportFName = '/home/geve/Dropbox/EnKF_advection/Run/advect|VAR|.eps'
$!Export 
  ExportRegion = AllFrames

$!Pick AddAtPosition
  X = 0.515570342205
  Y = 0.904524714829
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!Pick Clear
$!ENDLOOP


