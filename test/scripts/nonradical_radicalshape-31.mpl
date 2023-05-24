read "fglm_build_matrixn.mpl":

char:= nextprime (2^30):
vars:= [x,y]:
F:= subs ({x=4578*x+789*y,y=1453*x+977*y},[x^2,x*y,y^2]):
str:= "nonradical_radicalshape-31":

main (F,vars,char,str);
