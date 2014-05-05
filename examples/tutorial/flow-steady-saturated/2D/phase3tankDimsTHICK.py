#
# This simple python script generates the 2D geometry (Amanzi "regions") that describes the 2D
# "Closed Tank" scenario proposed by Greg Flach in a series of Emails to the Amanzi tank working
# group around 4/16/2014.  All regions are described using Amanzi "Polygons".  Additional
# logical regions will be used to group these by material property, but that is done outside this
# script.  A similar process will be used to generate 3D geometries based on the "Rotated Polygon"
# primitive.
#
#  MSD, 4/22/2014
#
    
# A handy function to print each region
def printRegion(coords,name):
    print "    <ParameterList name=\""+name+"\">"
    print "      <ParameterList name=\"Region: Polygon\">"
    print "        <Parameter name=\"VerticesV1\" type=\"Array(double)\" value=\"{",",".join([str(a[0]) for a in coords]),"}\"/>"
    print "        <Parameter name=\"VerticesV2\" type=\"Array(double)\" value=\"{",",".join([str(a[1]) for a in coords]),"}\"/>"
    print "      </ParameterList>"
    print "    </ParameterList>"

# Key dimensions
PrimaryLinerInsideDiameter = 24.5
PrimaryLinerInsideHeight   = 8
PrimaryLinerThickness      = .125
PrimaryPadThickness        = .25

SecondaryLinerAnnulus      = 0.75
SecondaryLinerThickness    = .125
SecondaryLinerInsideHeight = 1.5
SecondaryPadThickness      = .25

ConcreteFloorThickness     = 0.9
ConcreteRoofThickness      = 0.9
ConcreteWallThickness      = 0.7

SoilAboveTank              = 10
SoilBelowTank              = 10
SoilBesideTank             = 10 # Unused

FastFlowThickness          = .25
TankWasteLayerThickness    = .125

# Horizontal keypoints, from the center out
xDomain = 50
xCenter = xDomain/2
xInsidePrimaryLiner    = xCenter                - PrimaryLinerInsideDiameter/2
xOutsidePrimaryLiner   = xInsidePrimaryLiner    - PrimaryLinerThickness
xInsideAnnulus         = xOutsidePrimaryLiner   - SecondaryLinerAnnulus
xInsideSecondaryLiner  = xInsideAnnulus         - FastFlowThickness
xOutsideSecondaryLiner = xInsideSecondaryLiner  - SecondaryLinerThickness
xOutsideFF             = xOutsideSecondaryLiner - FastFlowThickness
xLeft                  = xOutsideSecondaryLiner - ConcreteWallThickness
xRight                 = xCenter - xLeft + xCenter

# Vertical keypoints, from the bottom up
yBottom = 0
yStart = yBottom + SoilBelowTank
yTopConcrete           = yStart                 + ConcreteFloorThickness
yTopSecondaryPad       = yTopConcrete           + SecondaryPadThickness
yTopSecondaryLiner     = yTopSecondaryPad       + SecondaryLinerThickness
yTopPrimaryPad         = yTopSecondaryLiner     + PrimaryPadThickness
yTopPrimaryLiner       = yTopPrimaryPad         + PrimaryLinerThickness
yTopWasteLayer         = yTopPrimaryLiner       + TankWasteLayerThickness
yTopSecondaryLinerLip  = yTopSecondaryLiner     + SecondaryLinerInsideHeight
yTopFF                 = yTopSecondaryLinerLip  + FastFlowThickness
yBottomInsidePrimaryLiner = yTopPrimaryLiner    + PrimaryLinerInsideHeight
yTopPrimaryLinerTOT    = yBottomInsidePrimaryLiner + PrimaryLinerThickness
yTop                   = yTopPrimaryLinerTOT    + ConcreteRoofThickness


# Exterior fast flow path, and its mirror image
FFLeft = []
FFLeft.append([xLeft,yTopConcrete])
FFLeft.append([xOutsideSecondaryLiner,yTopConcrete])
FFLeft.append([xOutsideSecondaryLiner,yTopSecondaryLinerLip])
FFLeft.append([xInsideSecondaryLiner,yTopSecondaryLinerLip])
FFLeft.append([xInsideSecondaryLiner,yTopPrimaryPad])
FFLeft.append([xInsideAnnulus,yTopPrimaryPad])
FFLeft.append([xInsideAnnulus,yTopFF])
FFLeft.append([xOutsideFF,yTopFF])
FFLeft.append([xOutsideFF,yTopConcrete+FastFlowThickness])
FFLeft.append([xLeft,yTopConcrete+FastFlowThickness])

FFRight = [[2*xCenter-a[0],a[1]] for a in FFLeft[::-1]]

# Interior fast flow path, and its mirror image
FFleft1 = []
FFleft1.append([xOutsidePrimaryLiner,yTopPrimaryPad])
FFleft1.append([xInsidePrimaryLiner,yTopPrimaryPad])
FFleft1.append([xInsidePrimaryLiner,yTopWasteLayer])
FFleft1.append([xOutsidePrimaryLiner,yTopWasteLayer])

FFright1 = [[2*xCenter-a[0],a[1]] for a in FFleft1[::-1]]


# Tank wall side, and its mirror image
WallLeft = []
WallLeft.append([xLeft,yTopConcrete+FastFlowThickness])
WallLeft.append([xOutsideFF,yTopConcrete+FastFlowThickness])
WallLeft.append([xOutsideFF,yTopFF])
WallLeft.append([xOutsideSecondaryLiner,yTopFF])
WallLeft.append([xOutsideSecondaryLiner,yTopPrimaryLinerTOT])
WallLeft.append([xLeft,yTopPrimaryLinerTOT])

WallRight = [[2*xCenter-a[0],a[1]] for a in WallLeft[::-1]]


# Tank top and bottom
WallTop = []
WallTop.append([xLeft,yTopPrimaryLinerTOT])
WallTop.append([xRight,yTopPrimaryLinerTOT])
WallTop.append([xRight,yTop])
WallTop.append([xLeft,yTop])

WallBottom = []
WallBottom.append([xLeft,yStart])
WallBottom.append([xRight,yStart])
WallBottom.append([xRight,yTopConcrete])
WallBottom.append([xLeft,yTopConcrete])

# Sand pads under primary and secondary liners
Pad1 = []
Pad1.append([xInsidePrimaryLiner,yTopSecondaryLiner])
Pad1.append([2*xCenter-xInsidePrimaryLiner,yTopSecondaryLiner])
Pad1.append([2*xCenter-xInsidePrimaryLiner,yTopSecondaryLiner+SecondaryPadThickness])
Pad1.append([xInsidePrimaryLiner,yTopSecondaryLiner+SecondaryPadThickness])

Pad2 = []
Pad2.append([xOutsideSecondaryLiner,yTopConcrete])
Pad2.append([2*xCenter-xOutsideSecondaryLiner,yTopConcrete])
Pad2.append([2*xCenter-xOutsideSecondaryLiner,yTopConcrete+PrimaryPadThickness])
Pad2.append([xOutsideSecondaryLiner,yTopConcrete+PrimaryPadThickness])


# Primary liner (bottom piece, and side/top piece
Liner1 = []
Liner1.append([xInsidePrimaryLiner,yTopPrimaryPad])
Liner1.append([2*xCenter - xInsidePrimaryLiner,yTopPrimaryPad])
Liner1.append([2*xCenter - xInsidePrimaryLiner,yTopPrimaryPad+PrimaryLinerThickness])
Liner1.append([xInsidePrimaryLiner,yTopPrimaryPad+PrimaryLinerThickness])

Liner1a = []
Liner1a.append([xOutsidePrimaryLiner,yTopWasteLayer])
Liner1a.append([xInsidePrimaryLiner,yTopWasteLayer])
Liner1a.append([xInsidePrimaryLiner,yBottomInsidePrimaryLiner])
Liner1a.append([2*xCenter-xInsidePrimaryLiner,yBottomInsidePrimaryLiner])
Liner1a.append([2*xCenter-xInsidePrimaryLiner,yTopWasteLayer])
Liner1a.append([2*xCenter-xOutsidePrimaryLiner,yTopWasteLayer])
Liner1a.append([2*xCenter-xOutsidePrimaryLiner,yTopPrimaryLinerTOT])
Liner1a.append([xOutsidePrimaryLiner,yTopPrimaryLinerTOT])

# Secondary liner
Liner2 = []
Liner2.append([xOutsideSecondaryLiner,yTopSecondaryPad])
Liner2.append([2*xCenter-xOutsideSecondaryLiner,yTopSecondaryPad])
Liner2.append([2*xCenter-xOutsideSecondaryLiner,yTopSecondaryLinerLip])
Liner2.append([2*xCenter-xInsideSecondaryLiner,yTopSecondaryLinerLip])
Liner2.append([2*xCenter-xInsideSecondaryLiner,yTopSecondaryLiner])
Liner2.append([xInsideSecondaryLiner,yTopSecondaryLiner])
Liner2.append([xInsideSecondaryLiner,yTopSecondaryLinerLip])
Liner2.append([xOutsideSecondaryLiner,yTopSecondaryLinerLip])

# Waste layer inside primary liner, and annular section (and its mirror image)
Waste1 = []
Waste1.append([xInsidePrimaryLiner,yTopPrimaryPad+PrimaryLinerThickness])
Waste1.append([2*xCenter-xInsidePrimaryLiner,yTopPrimaryPad+PrimaryLinerThickness])
Waste1.append([2*xCenter-xInsidePrimaryLiner,yTopWasteLayer])
Waste1.append([xInsidePrimaryLiner,yTopWasteLayer])

WasteL = []
WasteL.append([xInsideSecondaryLiner,yTopSecondaryLiner])
WasteL.append([xInsidePrimaryLiner,yTopSecondaryLiner])
WasteL.append([xInsidePrimaryLiner,yTopPrimaryPad])
WasteL.append([xInsideSecondaryLiner,yTopPrimaryPad])

WasteR = [[2*xCenter-a[0],a[1]] for a in WasteL[::-1]]

# A block encompassing the entire tank assembly
Tank = []
Tank.append([xLeft,yStart])
Tank.append([xRight,yStart])
Tank.append([xRight,yTop])
Tank.append([xLeft,yTop])

# Write out each region
printRegion(WallTop,"Wall Top")
printRegion(WallBottom,"Wall Bottom")
printRegion(WallRight,"Wall Right")
printRegion(WallLeft,"Wall Left")
printRegion(FFLeft,"FF Left")
printRegion(FFRight,"FF Right")
printRegion(Pad1,"Primary Sand Pad")
printRegion(Pad2,"Secondary Sand Pad")
printRegion(Liner1,"Primary Liner Floor")
printRegion(Liner1a,"Primary Liner")
printRegion(Liner2,"Secondary Liner")
printRegion(Waste1,"Waste Floor")
printRegion(WasteL,"Waste Left")
printRegion(WasteR,"Waste Right")
printRegion(FFleft1,"FF Left 1")
printRegion(FFright1,"FF Right 1")
printRegion(Tank,"Tank Block")
