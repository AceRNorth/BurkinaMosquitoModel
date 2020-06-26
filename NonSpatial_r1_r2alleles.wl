(* ::Package:: *)


(*homing distortion - em amd ef are homing rate 
parameters in males and females respectively*)
emplus=(1+em)/2;
efplus=(1+ef)/2;

(*sex-distortion - b=0 is default of no distortion*)
bplus=(1+b)/2;
bminus=(1-b)/2;

(*r2 generation from WH zygotes*)
\[Gamma]mdash=\[Gamma]m/2;
\[Gamma]fdash=\[Gamma]f/2;

(*wildtype generation from WH zygoes*)
\[Theta]m=(1-em-\[Gamma]m)/2;
\[Theta]f=(1-ef-\[Gamma]f)/2;
(*Fitnesses*)
\[Omega]whm=(1-\[Rho]m)(1-\[Epsilon]);
\[Omega]whf=(1-\[Rho]f)(1-\[Epsilon]);
\[Omega]whb=(1-\[Rho]f)(1-\[Rho]m)(1-\[Epsilon]);
\[Omega]wwm=\[Omega]wr1m=\[Omega]r1r1m=\[Omega]wr2m=\[Omega]r2r2m=\[Omega]r1r2m=(1-\[Rho]m);
\[Omega]wwf=\[Omega]wr1f=\[Omega]r1r1f=\[Omega]wr2f=\[Omega]r2r2f=\[Omega]r1r2f=(1-\[Rho]f);
\[Omega]wwb=\[Omega]wr1b=\[Omega]r1r1b=\[Omega]wr2b=\[Omega]r2r2b=\[Omega]r1r2b=(1-\[Rho]m)(1-\[Rho]f);
\[Omega]hr1m=(1-\[Rho]m)(1-\[Epsilon]);
\[Omega]hr1f=(1-\[Rho]f)(1-\[Epsilon]);
\[Omega]hr1b=(1-\[Rho]f)(1-\[Rho]m)(1-\[Epsilon]);

(*Male gamete frequencies*)
mwX=Mww/2+Mwr1/4+Mwr2/4;
mwXdash=bminus \[Theta]m Mwh;
mwY=Mww/2+Mwr1/4+Mwr2/4;
mwYdash=bplus \[Theta]m Mwh;
mhXdash=bminus emplus Mwh+bminus Mhh+bminus Mhr1/2+bminus Mhr2/2;
mhYdash=bplus emplus Mwh+bplus Mhh+bplus Mhr1/2+bplus Mhr2/2;
mr1X=Mwr1/4+Mr1r1/2+Mr1r2/4;
mr1Xdash=bminus Mhr1/2;
mr2X=Mwr2/4+Mr2r2/2+Mr1r2/4;
mr2Xdash=bminus \[Gamma]mdash Mwh+bminus Mhr2/2;
mr1Y=Mwr1/4+Mr1r1/2+Mr1r2/4;
mr1Ydash=bplus Mhr1/2;
mr2Y=Mwr2/4+Mr2r2/2+Mr1r2/4;
mr2Ydash=bplus \[Gamma]mdash Mwh+bplus Mhr2/2;

(*Female gamete frequencies. \[Omega]f and \[Omega]m are the fitnesses of 
w/d females where the d is maternal and paternal respectively*)
fw=(Fww+\[Omega]wwm FwwM+\[Omega]wwf FwwF+\[Omega]wwb FwwB)+(Fwr1+\[Omega]wr1m Fwr1M+\[Omega]wr1f Fwr1F+\[Omega]wr1b Fwr1B)/2+(Fwr2+\[Omega]wr2m Fwr2M+\[Omega]wr2f Fwr2F+\[Omega]wr2b Fwr2B)/2;

fwdash=\[Theta]f (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB);
fhdash=efplus (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB)+ (\[Omega]hr1f Fhr1F+\[Omega]hr1m Fhr1M+\[Omega]hr1b Fhr1B)/2;
fr1=(Fr1r1+\[Omega]r1r1m Fr1r1M+\[Omega]r1r1f Fr1r1F+\[Omega]r1r1b Fr1r1B)+
(Fwr1+\[Omega]wr1m Fwr1M+\[Omega]wr1f Fwr1F+\[Omega]wr1b Fwr1B)/2+
(Fr1r2+\[Omega]r1r2m Fr1r2M+\[Omega]r1r2f Fr1r2F+\[Omega]r1r2b Fr1r2B)/2;
fr1dash= (\[Omega]hr1f Fhr1F+\[Omega]hr1m Fhr1M+\[Omega]hr1b Fhr1B)/2;
fr2=(Fwr2+\[Omega]wr2m Fwr2M+\[Omega]wr2f Fwr2F+\[Omega]wr2b Fwr2B)/2+
(Fr1r2+\[Omega]r1r2m Fr1r2M+\[Omega]r1r2f Fr1r2F+\[Omega]r1r2b Fr1r2B)/2;
fr2dash=\[Gamma]fdash (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB);

(*next generation frequencies*)
nFww=mwX fw;
nFwwM=mwXdash fw;
nFwwF=mwX fwdash;
nFwwB=mwXdash fwdash;
nFwhM=mhXdash fw;
nFwhF=mwX fhdash;
nFwhB=mhXdash fwdash+mwXdash fhdash;
nFhh=mhXdash fhdash;
nFwr1=mr1X fw+mwX fr1;
nFwr1M=mr1Xdash fw+mwXdash fr1;
nFwr1F=mr1X fwdash+mwX fr1dash;
nFwr1B=mr1Xdash fwdash+mwXdash fr1dash;
nFwr2=mr2X fw+mwX fr2;
nFwr2M=mr2Xdash fw+mwXdash fr2;
nFwr2F=mr2X fwdash+mwX fr2dash;
nFwr2B=mr2Xdash fwdash+mwXdash fr2dash;
nFr1r1=mr1X fr1;
nFr1r1M=mr1Xdash fr1;
nFr1r1F=mr1X fr1dash;
nFr1r1B=mr1Xdash fr1dash;
nFr2r2=mr2X fr2+mr2Xdash fr2+mr2X fr2dash+mr2Xdash fr2dash;
nFhr1M=mhXdash fr1;
nFhr1F=mr1X fhdash;
nFhr1B=mhXdash fr1dash+mr1Xdash fhdash;
nFhr2=mhXdash fr2+mr2X fhdash+mhXdash fr2dash+mr2Xdash fhdash;
nFr1r2=mr1X fr2+mr1X fr2;
nFr1r2M=mr2Xdash fr1+mr1Xdash fr2;
nFr1r2F=mr2X fr1dash+mr1X fr2dash;
nFr1r2B=mr2Xdash fr1dash+mr1Xdash fr2dash;
nMww=(mwY+mwYdash)(fw+fwdash);
nMwh=mhYdash (fw+fwdash)+(mwY+mwYdash) fhdash;
nMhh=mhYdash fhdash;
nMwr1=(mr1Y+mr1Ydash) (fw+fwdash)+(mwY+mwYdash) (fr1+fr1dash);
nMr1r1=(mr1Y+mr1Ydash)(fr1+fr1dash);
nMhr1=mhYdash (fr1+fr1dash)+(mr1Y+mr1Ydash) fhdash;
nMwr2=(mr2Y+mr2Ydash) (fw+fwdash)+(mwY+mwYdash) (fr2+fr2dash);
nMr2r2=(mr2Y+mr2Ydash)(fr2+fr2dash);
nMhr2=mhYdash (fr2+fr2dash)+(mr2Y+mr2Ydash) fhdash;
nMr1r2=(mr2Y+mr2Ydash) (fr1+fr1dash)+(mr1Y+mr1Ydash) (fr2+fr2dash);

tot=Total[{nFww,nFwwM,nFwwF,nFwwB,nFwhM,nFwhF,nFwhB,nFhh,
nFwr1,nFwr1M,nFwr1F,nFwr1B,nFr1r1,nFr1r1M,nFr1r1F,nFr1r1B,nFhr1M,nFhr1F,nFhr1B,nFwr2,nFwr2M,
nFwr2F,nFwr2B,nFr2r2,nFhr2,nFr1r2,nFr1r2M,nFr1r2F,nFr1r2B,nMww,nMwh,nMhh,nMwr1,nMr1r1,nMhr1,nMwr2,nMr2r2,nMhr2,nMr1r2}];

next={nFww,nFwwM,nFwwF,nFwwB,nFwhM,nFwhF,nFwhB,nFhh,nFwr1,nFwr1M,nFwr1F,nFwr1B,
nFr1r1,nFr1r1M,nFr1r1F,nFr1r1B,nFhr1M,nFhr1F,nFhr1B,nFwr2,nFwr2M,nFwr2F,nFwr2B,
nFr2r2,nFhr2,nFr1r2,nFr1r2M,nFr1r2F,nFr1r2B,nMww,nMwh,nMhh,nMwr1,nMr1r1,nMhr1,nMwr2,nMr2r2,nMhr2,nMr1r2}/tot;



Runit[pars_,initial_,generations_]:=
Block[{list,lll,allelelist},
list=
NestList[next/.pars/.{Fww->#[[1]],FwwM->#[[2]],FwwF->#[[3]],FwwB->#[[4]],FwhM->#[[5]],FwhF->#[[6]],FwhB->#[[7]],Fhh->#[[8]],
Fwr1->#[[9]],Fwr1M->#[[10]],Fwr1F->#[[11]],Fwr1B->#[[12]],Fr1r1->#[[13]],Fr1r1M->#[[14]],
Fr1r1F->#[[15]],Fr1r1B->#[[16]],Fhr1M->#[[17]],Fhr1F->#[[18]],Fhr1B->#[[19]],Fwr2->#[[20]],
Fwr2M->#[[21]],Fwr2F->#[[22]],Fwr2B->#[[23]],Fr2r2->#[[24]],Fhr2->#[[25]],Fr1r2->#[[26]],
Fr1r2M->#[[27]],Fr1r2F->#[[28]],Fr1r2B->#[[29]],Mww->#[[30]],Mwh->#[[31]],Mhh->#[[32]],
Mwr1->#[[33]],Mr1r1->#[[34]],Mhr1->#[[35]],Mwr2->#[[36]],Mr2r2->#[[37]],Mhr2->#[[38]],Mr1r2->#[[39]]}&,initial,generations];
list];



