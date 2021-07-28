cd /Users/matteo/Coding/protein_modelling/ITPR3/ITPR3_analysis/structures

color white, not chain A and element C
color 0x40e0d0, chain A and element C
select itpr3, chain A and resi 196+2506+2524+615
select itpr1, chain A and resi 2005+2473
color 0xfa8072, itpr3 and element C
color 0xffd700, itpr1 and element C
show sticks, itpr3 or itpr1
set cartoon_transparency, 0.5, polymer and not chain A
show spheres, resn DUM
color grey60, resn DUM
backless
set_view (\
     0.622179210,   -0.050072074,    0.781272531,\
    -0.782290339,   -0.001253415,    0.622912109,\
    -0.030210685,   -0.998742998,   -0.039951675,\
     0.000009000,    0.000697648, -571.055419922,\
   212.546997070,  210.331649780,  205.474395752,\
   412.981597900,  729.135253906,  -20.000000000 )
   
set_color turquoise, [64,224,208]
set_color teal, [0, 128, 128]
set_color verdigris, [67, 179, 174]
set_color salmon, [255, 153, 153]
set_color aquamarine, [127, 255, 212]
set_color azure, [240, 255, 255]
set_color robin, [0, 204, 204]
set_color celeste, [178, 255, 255]
set_color lcyan, [224, 255, 255]
set_color lavender, [230, 230, 250]


## R colors
set_color r0, [248, 118, 109]
set_color r1, [196, 154, 0]
set_color r2, [83, 180, 0]
set_color r3, [0, 192, 148]
set_color r4, [0, 182, 235]
set_color r5, [165, 138, 255]
set_color r6, [251, 97, 215]

select cytoplasmic, chain A and resi 1-2202
select transmembrane, chain A and resi 2206-2671

select coupling_barrel, chain A and resi 3-230
select ITP_binding_barrel, chain A and resi 233-433
select first_ITP_binding_bundle, chain A and resi 437-707
select second_ITP_binding_bundle, chain A and resi 1175-1334
select coupling_bundle, chain A and resi 1335-1546
select above_membrane, chain A and resi 1864-1974
select ion_channel, chain A and resi 2191-2537

color white, not chain A and element C
color grey80, chain A and element C
color r6, chain A and element C and coupling_barrel
color r3, chain A and element C and ITP_binding_barrel
color r4, chain A and element C and first_ITP_binding_bundle
color r5, chain A and element C and second_ITP_binding_bundle
color r2, chain A and element C and coupling_bundle
color r1, chain A and element C and above_membrane
color r0, chain A and element C and ion_channel
color r6, element C and resi 3-230 and chain D
color r2, element C and resi 1335-1546 and chain B

set sphere_transparency, 0.5, resn DUM

hide sticks, not chain E+F+G+A
hide spheres, (not chain E+F+G) and not resn DUM
hide cartoon, (not chain A) and (not (resi 1335-1546 and chain B)) and (not (resi 3-230 and chain D))
set_view (\
     0.511221349,    0.752922297,   -0.414444149,\
     0.836638391,   -0.325625837,    0.440453827,\
     0.196672603,   -0.571910143,   -0.796393275,\
    -0.003035948,   -0.000205636, -636.699340820,\
   238.736068726,  185.996185303,  235.230728149,\
  -1247.291870117, 2520.682128906,  -20.000000000 )

set cartoon_transparency, 0.3, polymer and not chain A
   
   
   
set_color r0, [248, 118, 109]
set_color r1, [124, 174, 0]
set_color r2, [0, 191, 196]
set_color r3, [199, 124, 255]
hide all
show cartoon
show sphere, resn DUM
set cartoon_transparency, 0, polymer
color r0, chain A
color r1, chain D
color r2, chain C
color r3, chain B
remove not (chain A+B+C+D or resn DUM)
show surface
set_view (\
     0.876986265,    0.159974396,   -0.453105658,\
     0.480366379,   -0.268665135,    0.834902465,\
     0.011827621,   -0.949856699,   -0.312463999,\
    -0.001747698,    0.000810385, -561.165954590,\
   208.144470215,  214.182235718,  218.605468750,\
   403.015747070,  719.169494629,  -20.000000000 )

##
cd /Users/matteo/Coding/protein_modelling/ITPR3/ITPR3_analysis/structures
set_color turquoise, [64,224,208]
set_color teal, [0, 128, 128]
set_color verdigris, [67, 179, 174]
set_color salmon, [255, 153, 153]
set_color aquamarine, [127, 255, 212]
set_color azure, [240, 255, 255]
set_color robin, [0, 204, 204]
set_color celeste, [178, 255, 255]
set_color lcyan, [224, 255, 255]
set_color lavender, [230, 230, 250]
load 6DRC.min.al.cif, wildtype
load variants/6DRC.A196T.pdb, mutant
align mutant and resi 195-197 and chain A, wildtype and resi 195-197 and chain A
color turquoise, wildtype and element C
color salmon, mutant and element C
color white, (not chain A) and element C
zoom resi 196 and chain A
show sticks, resi 219+38+26+196 and chain A
hide sticks, element H
set_view (\
     0.072623543,    0.970286906,    0.230805382,\
     0.995787799,   -0.083526760,    0.037810326,\
     0.055967245,    0.227087379,   -0.972264767,\
     0.000289530,   -0.000241246,  -55.646717072,\
   272.538208008,  237.100860596,  162.566513062,\
    40.421314240,   70.869842529,  -20.000000000 )

##
cd /Users/matteo/Coding/protein_modelling/ITPR3/ITPR3_analysis/structures
set_color turquoise, [64,224,208]
set_color teal, [0, 128, 128]
set_color verdigris, [67, 179, 174]
set_color salmon, [255, 153, 153]
set_color aquamarine, [127, 255, 212]
set_color azure, [240, 255, 255]
set_color robin, [0, 204, 204]
set_color celeste, [178, 255, 255]
set_color lcyan, [224, 255, 255]
set_color lavender, [230, 230, 250]
load 6DRC.min.al.cif, wildtype
load variants/6DRC.R2524H.pdb, mutant
align mutant and resi 2523-2525 and chain A, wildtype and resi 2523-2525 and chain A
color turquoise, wildtype and element C
color salmon, mutant and element C
color white, (not chain A) and element C
select variant_resi, resi 2524 and chain A
zoom variant_resi
## show sticks, byres (variant_resi expand 3)
show sticks, (resi 2524+2528 and chain A) or (resi 2522+2518 and chain D)
hide sticks, element H
    
    
##

cd /Users/matteo/Coding/protein_modelling/ITPR3/ITPR3_analysis/structures
set_color turquoise, [64,224,208]
set_color teal, [0, 128, 128]
set_color verdigris, [67, 179, 174]
set_color salmon, [255, 153, 153]
set_color aquamarine, [127, 255, 212]
set_color azure, [240, 255, 255]
set_color robin, [0, 204, 204]
set_color celeste, [178, 255, 255]
set_color lcyan, [224, 255, 255]
set_color lavender, [230, 230, 250]
load 6DRC.min.al.cif, wildtype
load variants/6DRC.I2506N.pdb, mutant
align mutant and resi 2505-2507 and chain A, wildtype and resi 2505-2507 and chain A
color turquoise, wildtype and element C
color salmon, mutant and element C
color white, (not chain A) and element C
select variant_resi, resi 2506 and chain A
select gate, resi 2517+2513+2472+2478
show sticks, gate and element C
zoom variant_resi
## show sticks, not element H and byres (variant_resi expand 3)
show sticks, (resi 2506+2510 and chain A)
hide sticks, element H