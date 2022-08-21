# Bone-Analysis

The geometric morphometric (GM) protocol captures the three-dimensional shape of the femur using 12 landmarks (Table 3, Figure 1) and four curves represented by 31 semi-landmarks. Eleven of the twelve landmarks and two curves (the 'femoral head' and 'intertrochanteric crest' curves) were positioned manually on the 3D surface model using EVAN Toolbox (ET)307. One landmark ('centroid of the femoral head') and the other two curves (the 'femoral neck' and 'subtrochanteric' curves) were created semi-automatically with a dedicated software we developed in MATLAB R2013a hereby provided.

The protocol was applied to left femora only. When the left femur was unavailable (e.g., fractured), the right femur was mirrored and used instead. Bilateral symmetry was examined to ensure the validity of using the contralateral femur.

Femoral neck:

Semi-automatic curve placed around the femoral neck at its midpoint. An anchor point is placed at landmark 8, and then two additional points are chosen on the anterior and posterior aspects of the femoral neck at its midpoint. The program calculates the area of all curves that pass through the anchor point and ±15° from the two selected points, and then chooses the curve with the smallest cross-sectional area.

Subtrochanteric:

Semi-automatic curve placed around the femoral shaft below the lesser trochanter.  An anchor point is placed at landmark 6, and then two additional landmarks are placed on the medial and lateral aspects of the femoral shaft at the same level of the anchor point. The program calculates the area of all curves that pass through the anchor point and ±15° from the two selected points, and then chooses the curve with the smallest cross-sectional area.

Instructions:

1. Run file 'CurveProgram30_Sphere", select proximal femur .stl
2. Select five points along on the femoral head (left click + enter) -> center of femoral head is provided (generate landmark in Evan Toolbox) 
3. Select 'one point' for the 'Find Curve' function of the program. Left click to select one point on the curve, then use left click to select two additional points. 
4. .Obj is generated -> use to generate curve in Evan Toolbox
