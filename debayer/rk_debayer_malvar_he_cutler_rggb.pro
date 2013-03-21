;+
; IDL IMPLEMENTATION FOR MALVAR HE CUTLER INTERPOLATION ALGORITHM 
; NAME:
;    rk_debayer_malvar_he_cutler_rggb
;
; PURPOSE:
;    This function takes an array with an R, G1, G2, B bayer pattern, where
;    0,0 is R and returns a combined array containing the debayered arrays R, G,
;    and B.
;
; REFERENCE:
;    High-quality linear interpolation for demosaicing of Bayer-patterned color images
;    By Henrique S. Malvar, Li-wei He, and Ross Cutler
;      http://research.microsoft.com/~rcutler/pub/Demosaicing_ICASSP04.pdf
;    RK_DEBAYER_MALVAR_HE_CUTLER_RGGB.PRO may be linked to U.S. Patent 7502505 by 
;    Henrique S. Malvar, Li-wei He, and Ross Cutler
;
; CATEGORY:
;   Interpolation algorithm
;
; CALLING SEQUENCE:
;   output_channels = rk_debayer_malvar_he_cutler_rggb(input_image_array)
;
; INPUTS:
;   input_image_array : RGGB Bayer Pattern Image in Array Format
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;   output_channels = [[[R]], [[G]], [[B]]]
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;    Currently only implemented for RGGB Bayer Pattern.
;
; PROCEDURE:
;
;
; TODO:
;
; EXAMPLE:
;   IDL> help, a
;   A               BYTE      = Array[1344, 1200]
;   IDL> b = rk_debayer_malvar_he_cutler_rggb(a)
;   IDL> help, b
;   B               INT       = Array[1344, 1200, 3]
;
; MODIFICATION HISTORY:
;    Written by: Rohan Kadam, 03/20/2013.
;
; NOTE:
;    Copyright (C) 2013 by rk
;    This library is free software; you can redistribute it and/or
;    modify.
;    This library is distributed in the hope that it will be useful for
;    science and education but WITHOUT ANY WARRANTY that it will always
;    fit the requirement.
;    See the REFERENCE above for other possible restrictions.
; 
;-
FUNCTION rk_debayer_malvar_he_cutler_rggb, IMAGE_ARRAY
  ;get input image array size
  image_array_size = size(image_array)

  ;input image width and height
  image_height = image_array_size[1]
  image_width = image_array_size[2]

  ;create a bordered input image to handle boundary pixel conditions
  b_input_image = replicate(0, image_height+4, image_width+4)

  ;copy input image array to bordered image
  b_input_image[2,2] = image_array

  ;create output matrices of the same size as the input matrix
  R_image_matrix = indgen(image_height,image_width)
  G_image_matrix = indgen(image_height,image_width)
  B_image_matrix = indgen(image_height,image_width)

  ;Initialize all the values to 0
  R_image_matrix[*,*] = 0
  G_image_matrix[*,*] = 0
  B_image_matrix[*,*] = 0

  ;loop through the bordered image matrix starting at 2,2 as the input image starts at 2,2
  for x = 2, image_height + 1 do begin
      for y = 2, image_width + 1 do begin
          if(x mod 2 EQ 0) then begin
              if(y mod 2 EQ 0) then begin
                  ;The pixel at x,y is a Red pixel
                  R_image_matrix[x-2,y-2] = b_input_image[x,y]
                  ;Calculate green pixel at x,y
                  G_image_matrix[x-2,y-2] =(                                   - b_input_image[x-2,y] $
                                                                              + 2*b_input_image[x-1,y] $
                                - b_input_image[x,y-2] + 2*b_input_image[x,y-1] + 4*b_input_image[x,y]   + 2*b_input_image[x,y+1] - b_input_image[x,y+2] $
                                                                              + 2*b_input_image[x+1,y] $
                                                                              - b_input_image[x+2,y] $
                                           )/8
                  ;Calculate blue value at x,y
                  B_image_matrix[x-2,y-2] = (                                - (3*b_input_image[x-2,y])/2 $
                                                  + 2*b_input_image[x-1,y-1] +        0                    + 2*b_input_image[x-1,y+1] $
                     - (3*b_input_image[x,y-2])/2 +          0               + 6*b_input_image[x,y]        +            0               - (3*b_input_image[x,y+2])/2 $
                                                  + 2*b_input_image[x+1,y-1] +        0                    + 2*b_input_image[x+1,y+1] $
                                                                             - (3*b_input_image[x+2,y])/2 $
                                           )/8

              endif else begin
                  ;green pixel in red row and blue column
                  G_image_matrix[x-2,y-2] = b_input_image[x,y]

                  ;Calculate Red value at x,y
                  R_image_matrix[x-2,y-2] = (                                 (1*b_input_image[x-2,y])/2 $
                                                  -   b_input_image[x-1,y-1] +        0                    - b_input_image[x-1,y+1] $
                           - b_input_image[x,y-2] + 4*b_input_image[x,y-1]   + 5*b_input_image[x,y]        + 4*b_input_image[x,y+1]    - b_input_image[x,y+2] $
                                                  - b_input_image[x+1,y-1]   +        0                    - b_input_image[x+1,y+1] $
                                                                             + (b_input_image[x+2,y])/2 $
                                           )/8
                  ;Calculate Green Value at x,y
                  B_image_matrix[x-2,y-2] =(                                 - b_input_image[x-2,y] $
                                                  -   b_input_image[x-1,y-1] +        0                    - b_input_image[x-1,y+1] $
                         + b_input_image[x,y-2]/2 + 4*b_input_image[x,y-1]   + 5*b_input_image[x,y]        + 4*b_input_image[x,y+1]    + b_input_image[x,y+2]/2 $
                                                  - b_input_image[x+1,y-1]   +        0                    - b_input_image[x+1,y+1] $
                                                                             -  b_input_image[x+2,y] $
                                           )/8
              endelse
          endif else begin
              if(y mod 2 EQ 1) then begin
                  ;blue pixel at x,y
                  B_image_matrix[x-2,y-2] = b_input_image[x,y]

                  ;Calculate red value at x,y
                  R_image_matrix[x-2,y-2] = (                                - (3*b_input_image[x-2,y])/2 $
                                                  + 2*b_input_image[x-1,y-1] +        0                    + 2*b_input_image[x-1,y+1] $
                     - (3*b_input_image[x,y-2])/2 +          0               + 6*b_input_image[x,y]        +            0               - (3*b_input_image[x,y+2])/2 $
                                                  + 2*b_input_image[x+1,y-1] +        0                    + 2*b_input_image[x+1,y+1] $
                                                                             - (3*b_input_image[x+2,y])/2 $
                                           )/8

                  ;Calculate green value at x,y
                  G_image_matrix[x-2,y-2] = (                                   - b_input_image[x-2,y] $
                                                                              + 2*b_input_image[x-1,y] $
                              - b_input_image[x,y-2] + 2*b_input_image[x,y-1] + 4*b_input_image[x,y]       + 2*b_input_image[x,y+1] - b_input_image[x,y+2] $
                                                                              + 2*b_input_image[x+1,y] $
                                                                              - b_input_image[x+2,y] $
                                           )/8
              endif else begin
                  ;green pixel in blue row and red column
                  G_image_matrix[x-2,y-2] = b_input_image[x,y]
                  ;Calculate red value at x,y
                  R_image_matrix[x-2,y-2] = (                                 - b_input_image[x-2,y] $
                                                    - b_input_image[x-1,y-1]  +        0                 - b_input_image[x-1,y+1] $
                           + b_input_image[x,y-2]/2 + 4*b_input_image[x,y-1]  + 5*b_input_image[x,y]     + 4*b_input_image[x,y+1]    + b_input_image[x,y+2]/2 $
                                                    - b_input_image[x+1,y-1]  +        0                 - b_input_image[x+1,y+1] $
                                                                              -  b_input_image[x+2,y] $
                                           )/8
                  ;Calculate green value at x,y
                  B_image_matrix[x-2,y-2] = (                                  (1*b_input_image[x-2,y])/2 $
                                                  -   b_input_image[x-1,y-1] +        0                    - b_input_image[x-1,y+1] $
                           - b_input_image[x,y-2] + 4*b_input_image[x,y-1]   + 5*b_input_image[x,y]        + 4*b_input_image[x,y+1]    - b_input_image[x,y+2] $
                                                  - b_input_image[x+1,y-1]   +        0                    - b_input_image[x+1,y+1] $
                                                                             + (b_input_image[x+2,y])/2 $
                                           )/8
              endelse
          endelse
      endfor
  endfor
  return,([[[R_image_matrix]], [[G_image_matrix]],[[B_image_matrix]]])
END
