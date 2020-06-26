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
\[Omega]wwm=\[Omega]wrm=(1-\[Rho]m);
\[Omega]wwf=\[Omega]wrf=(1-\[Rho]f);
\[Omega]wwb=\[Omega]wrb=(1-\[Rho]m)(1-\[Rho]f);

(*Male gamete frequencies; dash signifies that the parent carrying a drive allele*)
mwX=Mww/2+Mwr/4;
mwXdash=bminus \[Theta]m Mwh;
mwY=Mww/2+Mwr/4;
mwYdash=bplus \[Theta]m Mwh;
mhXdash=bminus emplus Mwh+bminus Mhh+bminus Mhr/2;
mhYdash=bplus emplus Mwh+bplus Mhh+bplus Mhr/2;
mrX=Mwr/4+Mrr/2;
mrXdash=bminus \[Gamma]mdash Mwh+bminus Mhr/2;
mrY=Mwr/4+Mrr/2;
mrYdash=bplus \[Gamma]mdash Mwh+bplus Mhr/2;

(*Female gamete frequencies. \[Omega]f and \[Omega]m are the fitnesses of 
w/d females where the d is maternal and paternal respectively*)
fw=(Fww+\[Omega]wwm FwwM+\[Omega]wwf FwwF+\[Omega]wwb FwwB)+(Fwr+\[Omega]wrm FwrM+\[Omega]wrf FwrF+\[Omega]wrb FwrB)/2;
fwdash=\[Theta]f (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB);
fhdash=efplus (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB);
fr=(Fwr+\[Omega]wrm FwrM+\[Omega]wrf FwrF+\[Omega]wrb FwrB)/2;
frdash=\[Gamma]fdash (\[Omega]whf FwhF+\[Omega]whm FwhM+\[Omega]whb FwhB);

(*next generation frequencies*)
nFww=mwX fw;
nFwwM=mwXdash fw;
nFwwF=mwX fwdash;
nFwwB=mwXdash fwdash;
nFwhM=mhXdash fw;
nFwhF=mwX fhdash;
nFwhB=mhXdash fwdash+mwXdash fhdash;
nFhh=mhXdash fhdash;
nFwr=mrX fw+mwX fr;
nFwrM=mrXdash fw+mwXdash fr;
nFwrF=mrX fwdash+mwX frdash;
nFwrB=mrXdash fwdash+mwXdash frdash;
nFrr=mrX fr+mrXdash fr+mrX frdash+mrXdash frdash;
nFhr=mhXdash fr+mrX fhdash+mhXdash frdash+mrXdash fhdash;
nMww=(mwY+mwYdash)(fw+fwdash);
nMwh=mhYdash (fw+fwdash)+(mwY+mwYdash) fhdash;
nMhh=mhYdash fhdash;
nMwr=(mrY+mrYdash) (fw+fwdash)+(mwY+mwYdash) (fr+frdash);
nMrr=(mrY+mrYdash)(fr+frdash);
nMhr=mhYdash (fr+frdash)+(mrY+mrYdash) fhdash;

(*Normalise and simplify vector*)

tot=Simplify[Total[{nFww,nFwwM,nFwwF,nFwwB,nFwhM,nFwhF,nFwhB,nFhh,nFwr,nFwrM,nFwrF,nFwrB,nFrr,nFhr,nMww,nMwh,nMhh,nMwr,nMrr,nMhr}]]
next=Simplify[{nFww,nFwwM,nFwwF,nFwwB,nFwhM,nFwhF,nFwhB,nFhh,nFwr,nFwrM,nFwrF,nFwrB,nFrr,nFhr,nMww,nMwh,nMhh,nMwr,nMrr,nMhr}/tot];

(*Append with genetic load*)
AppendTo[next,1-2(Fwr+FwrF+FwrM+Fww+FwwB+FwwF+FwwM+(FwhB+FwhF+FwhM) \[Epsilon]-
(FwrF+FwwB+FwwF+(FwhB+FwhF) \[Epsilon]) \[Rho]f+FwrB (-1+\[Rho]f) (-1+\[Rho]m)-(FwrM+FwwB+FwwM+(FwhB+FwhM) \[Epsilon]-(FwwB+FwhB \[Epsilon]) \[Rho]f) \[Rho]m)];

(*Run model for given number of generations; pars needs to be declared as a substitution list, e.g. {em\[Rule]0.95,ef\[Rule]0.95,b->0,\[Omega]f->0.65,\[Omega]m->0.217,\[Gamma]m\[Rule]0.025,\[Gamma]f\[Rule]0.025,\[Rho]m\[Rule]0.2,\[Rho]f\[Rule]0.4,\[Epsilon]\[Rule]0.3}*)
Runit[pars_,initial_,generations_]:=
Block[{list,lll,allelelist},
list=
NestList[next/.pars/.{Fww->#[[1]],FwwM->#[[2]],FwwF->#[[3]],FwwB->#[[4]],FwhM->#[[5]],
FwhF->#[[6]],FwhB->#[[7]],Fhh->#[[8]],Fwr->#[[9]],FwrM->#[[10]],FwrF->#[[11]],
FwrB->#[[12]],Frr->#[[13]],Fhr->#[[14]],Mww->#[[15]],Mwh->#[[16]],Mhh->#[[17]],
Mwr->#[[18]],Mrr->#[[19]],Mhr->#[[20]],tt->#[[21]]}&,initial,generations];
list];
(*Run model until equilibrium is reach (change in genetic proportions is below a tolerance threshold*)

RunitEq[pars_,initial_,threshold_,maxgen_]:=Block[{list,lll,vec,vecnew,gen,dif},
vec=initial;
vecnew=vec;
gen=0;
dif=100;
While[gen<maxgen && dif>threshold,
vecnew=next/.pars/.{Fww->vec[[1]],FwwM->vec[[2]],FwwF->vec[[3]],FwwB->vec[[4]],
FwhM->vec[[5]],FwhF->vec[[6]],FwhB->vec[[7]],Fhh->vec[[8]],Fwr->vec[[9]],
FwrM->vec[[10]],FwrF->vec[[11]],FwrB->vec[[12]],Frr->vec[[13]],
Fhr->vec[[14]],Mww->vec[[15]],Mwh->vec[[16]],Mhh->vec[[17]],
Mwr->vec[[18]],Mrr->vec[[19]],Mhr->vec[[20]],tt->vec[[21]]};
dif=Max[Abs[vec-vecnew]];
vec=Chop[vecnew];
gen++];
{gen,vec}];

