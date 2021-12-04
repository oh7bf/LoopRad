#!/usr/bin/wolfram -script
(* Calculate rectangular loop radiation fields. *)
(* Free space loop on xy-plane with z=0. *)
(* Loop width 2a and height 2b. *)
(* Loop center at x = 0, y = 0 and z = 0. *)
(* Constant RF current Irf[A] along circumference. *)
(* Integration path clockwise. *)
(*                                          *)
(*         (-a, +b, 0)  Irf   (+a, +b, 0)   *)
(*               -------->---------         *)
(*              |                  |        *)
(* ^ y          |                  |        *)
(* |            |                  |        *)
(* |--> x       |        C         |        *)
(*               --------||--------         *)
(*         (-a, -b, 0)        (+a, -b, 0)   *)           
(*                                          *)
(* Edit frequency f[Hz], loop current Irf[A], width 2b[m] and height 2a[m] *)
(* below for other geometries. N.B. in real loop antenna Irf varies along *)
(* loop circumference, but constant Irf is easier to calculate. *)

(* References:*)  
(* Kai Siwiak, KE4PT and Richard Quick, W4RQ, Small Gap-resonated HF Loop
Antenna Fed by a Secondary Loop,  QEX July/August 2018 *)
(* Warren L. Stutzman, Gary A. Thiele, Antenna Theory and Design, 1998 John Wiley & Sons. *)

(* Sat Dec  4 13:16:49 CST 2021 *)
(* Jaakko Koivuniemi, AC1BB *)

(* accuracy goal in numerical integration *)
aG = 5

(* print caculation results *)
printCalcs = False

(* calculate H vector field *)
calcHfield = True 

(* calculate E vector field *)
calcEfield = False

(* calculate S vector field *)
calcSfield = False

(* Vector potential *)
Print["-- Vector potentials"]
Ax[\[Beta]_, \[Mu]r_, Irf_, a_, b_, x_, y_, z_] := \[Mu]r \[Mu]0 Irf  NIntegrate[ Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y - b)^2 + z^2]]/(4 Pi Sqrt[(x - xl)^2 + (y - b)^2 + z^2]), {xl, -a, +a}, AccuracyGoal->aG ] + \[Mu]r \[Mu]0 Irf  NIntegrate[ Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y + b)^2 + z^2]]/(4 Pi Sqrt[(x - xl)^2 + (y + b)^2 + z^2]), {xl, +a, -a}, AccuracyGoal->aG ]
If[printCalcs, Print[ ?Ax ], Nothing ]

Ay[\[Beta]_, \[Mu]r_, Irf_, a_, b_, x_, y_, z_] := \[Mu]r \[Mu]0 Irf  NIntegrate[ Exp[-I \[Beta] Sqrt[(x - a)^2 + (y - yl)^2 + z^2]]/(4 Pi Sqrt[(x - a)^2 + (y - yl)^2 + z^2]), {yl, +b, -b}, AccuracyGoal->aG ] + \[Mu]r \[Mu]0 Irf  NIntegrate[ Exp[-I \[Beta] Sqrt[(x + a)^2 + (y - yl)^2 + z^2]]/(4 Pi Sqrt[(x + a)^2 + (y - yl)^2 + z^2]), {yl, -b, +b}, AccuracyGoal->aG ]
If[printCalcs, Print[ ?Ay ], Nothing ]

(* Magnetic field *)
Print["-- Magnetic field"]

(* Ax components *)
curlAx1 = Curl[{Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y - b)^2 + z^2]]/ Sqrt[(x - xl)^2 + (y - b)^2 + z^2], 0, 0}, {x, y, z}] 
Print[ "-- Curl Ax1 / constants" ]
If[printCalcs, Print[ curlAx1 ], Nothing ]

curlAx2 = Curl[{Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y + b)^2 + z^2]]/ Sqrt[(x - xl)^2 + (y + b)^2 + z^2], 0, 0}, {x, y, z}]
Print[ "-- Curl Ax2 / constants" ]
If[printCalcs, Print[ curlAx2 ], Nothing ]

(* Ay components *)
curlAy1 = Curl[{0, Exp[-I \[Beta] Sqrt[(x - a)^2 + (y - yl)^2 + z^2]]/ Sqrt[(x - a)^2 + (y - yl)^2 + z^2], 0}, {x, y, z}]
Print[ "-- Curl Ay1 / constants" ]
If[printCalcs, Print[ curlAy1 ], Nothing ]

curlAy2 = Curl[{0, Exp[-I \[Beta] Sqrt[(x + a)^2 + (y - yl)^2 + z^2]]/ Sqrt[(x + a)^2 + (y - yl)^2 + z^2], 0}, {x, y, z}]
Print[ "-- Curl Ay2 / constants" ]
If[printCalcs, Print[ curlAy2 ], Nothing ]

(* Hx integration *)
Hx[\[Beta]_, Irf_, a_, b_, x_, y_, z_] := (Irf/(4 Pi)) (NIntegrate[( E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) z)/((-a + x)^2 + (y - yl)^2 + z^2)^(3/2) + ( I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) z \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2) , {yl, +b, -b}, AccuracyGoal->aG ] + NIntegrate[( E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) z)/((a + x)^2 + (y - yl)^2 + z^2)^(3/2) + ( I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) z \[Beta])/((a + x)^2 + (y - yl)^2 + z^2), {yl, -b, +b}, AccuracyGoal->aG ] )
If[printCalcs, Print[ ?Hx ], Nothing ]

(* Hy integration *)
Hy[\[Beta]_, Irf_, a_, b_, x_, y_, z_] := (Irf/(4 Pi)) (NIntegrate[-(( E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) z)/((x - xl)^2 + (-b + y)^2 + z^2)^(3/2)) - ( I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) z \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2), {xl, -a, +a}, AccuracyGoal->aG ] + NIntegrate[-(( E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) z)/((x - xl)^2 + (b + y)^2 + z^2)^(3/2)) - ( I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) z \[Beta])/((x - xl)^2 + (b + y)^2 + z^2), {xl, +a, -a}, AccuracyGoal-> aG ] )
If[printCalcs, Print[ ?Hy ], Nothing ]

(* Hz integration *)
Hz[\[Beta]_, Irf_, a_, b_, x_, y_, z_] := (Irf/(4 Pi)) (NIntegrate[( E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (-b + y))/((x - xl)^2 + (-b + y)^2 + z^2)^(3/2) + ( I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (-b + y) \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2), {xl, -a, +a}, AccuracyGoal->aG ] + NIntegrate[( E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (b + y))/((x - xl)^2 + (b + y)^2 + z^2)^(3/2) + ( I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (b + y) \[Beta])/((x - xl)^2 + (b + y)^2 + z^2), {xl, +a, -a}, AccuracyGoal->aG ] + NIntegrate[-(( E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (-a + x))/((-a + x)^2 + (y - yl)^2 + z^2)^(3/2)) - ( I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (-a + x) \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2), {yl, +b, -b}, AccuracyGoal->aG ] + NIntegrate[-(( E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (a + x))/((a + x)^2 + (y - yl)^2 + z^2)^(3/2)) - ( I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (a + x) \[Beta])/((a + x)^2 + (y - yl)^2 + z^2), {yl, -b, +b}, AccuracyGoal->aG ] )
If[printCalcs, Print[ ?Hz ], Nothing ]

(* Htotal total magnetic field *)
Htotal[\[Beta]_, Irf_, a_, b_, x_, y_, z_] := Sqrt[Abs[Hx[\[Beta], Irf, a, b, x, y, z]]^2 + Abs[Hy[\[Beta], Irf, a, b, x, y, z]]^2 + Abs[Hz[\[Beta], Irf, a, b, x, y, z]]^2]
If[printCalcs, Print[ ?Htotal ], Nothing ]

(* Electric field *)

(* Ax components *)
gradDivAx1 = Grad[Div[{Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y - b)^2 + z^2]]/ Sqrt[(x - xl)^2 + (y - b)^2 + z^2], 0, 0}, {x, y, z}], {x, y, z}] 
Print[ "-- Grad Div Ax1 / constants" ]
If[printCalcs, Print[ gradDivAx1 ], Nothing ]

gradDivAx2 = Grad[Div[{Exp[-I \[Beta] Sqrt[(x - xl)^2 + (y - b)^2 + z^2]]/ Sqrt[(x - xl)^2 + (y - b)^2 + z^2], 0, 0}, {x, y, z}], {x, y, z}] 
Print[ "-- Grad Div Ax2 / constants" ]
If[printCalcs, Print[ gradDivAx2 ], Nothing ]

(* Ay components *)
gradDivAy1 = Grad[Div[{0, Exp[-I \[Beta] Sqrt[(x - a)^2 + (y - yl)^2 + z^2]]/ Sqrt[(x - a)^2 + (y - yl)^2 + z^2], 0}, {x, y, z}], {x, y, z}] 
Print[ "-- Grad Div Ay1 / constants" ]
If[printCalcs, Print[ gradDivAy1 ], Nothing ]

gradDivAy2 = Grad[Div[{0, Exp[-I \[Beta] Sqrt[(x + a)^2 + (y - yl)^2 + z^2]]/ Sqrt[(x + a)^2 + (y - yl)^2 + z^2], 0}, {x, y, z}], {x, y, z}]
Print[ "-- Grad Div Ay2 / constants" ]
If[printCalcs, Print[ gradDivAy2 ], Nothing ]

(* Ex integration *)
Ex[\[Omega]_, \[Beta]_, \[Mu]r_, \[Epsilon]r_, Irf_, a_, b_, x_, y_, 
  z_] := -I \[Omega] Ax[\[Beta], \[Mu]r, Irf, a, b, x, y, 
    z] - (I Irf/(4 Pi \[Omega]  \[Epsilon]r \[Epsilon]0)) \
(NIntegrate[(
       3 E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2)/((x - xl)^2 + (-b + y)^2 + z^2)^(5/2) - 
       E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + 
         z^2] \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2)^(3/2) + (
       3 I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2 \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2)^2 - (
       I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + 
          z^2] \[Beta]) \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2) - (
       E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2 \[Beta]^2)/((x - xl)^2 + (-b + y)^2 + z^2)^(
       3/2), {xl, -a, +a}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2)/((x - xl)^2 + (b + y)^2 + z^2)^(5/2) - 
       E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + 
         z^2] \[Beta])/((x - xl)^2 + (b + y)^2 + z^2)^(3/2) + (
       3 I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2 \[Beta])/((x - xl)^2 + (b + y)^2 + z^2)^2 - (
       I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + 
          z^2] \[Beta]) \[Beta])/((x - xl)^2 + (b + y)^2 + z^2) - (
       E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
          xl)^2 \[Beta]^2)/((x - xl)^2 + (b + y)^2 + z^2)^(
       3/2), {xl, +a, -a}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (-a + 
          x) (y - yl))/((-a + x)^2 + (y - yl)^2 + z^2)^(5/2) + (
       3 I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (-a + 
          x) (y - yl) \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2)^2 - (
       E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (-a + 
          x) (y - yl) \[Beta]^2)/((-a + x)^2 + (y - yl)^2 + z^2)^(
       3/2), {yl, +b, -b}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (a + 
          x) (y - yl))/((a + x)^2 + (y - yl)^2 + z^2)^(5/2) + (
       3 I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (a + 
          x) (y - yl) \[Beta])/((a + x)^2 + (y - yl)^2 + z^2)^2 - (
       E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (a + x) (y - 
          yl) \[Beta]^2)/((a + x)^2 + (y - yl)^2 + z^2)^(
       3/2), {yl, -b, +b}, AccuracyGoal->aG ] ) 

Ey[\[Omega]_, \[Beta]_, \[Mu]r_, \[Epsilon]r_, Irf_, a_, b_, x_, y_, 
  z_] := -I \[Omega] Ay[\[Beta], \[Mu]r, Irf, a, b, x, y, 
    z] - (I Irf/(4 Pi \[Omega]  \[Epsilon]r \[Epsilon]0)) \
(NIntegrate[(
       3 E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl) (-b + y))/((x - xl)^2 + (-b + y)^2 + z^2)^(5/2) + (
       3 I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl) (-b + y) \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2)^2 - (
       E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
          xl) (-b + y) \[Beta]^2)/((x - xl)^2 + (-b + y)^2 + z^2)^(
       3/2), {xl, -a, +a}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
          xl) (b + y))/((x - xl)^2 + (b + y)^2 + z^2)^(5/2) + (
       3 I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
          xl) (b + y) \[Beta])/((x - xl)^2 + (b + y)^2 + z^2)^2 - (
       E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - xl) (b +
           y) \[Beta]^2)/((x - xl)^2 + (b + y)^2 + z^2)^(
       3/2), {xl, +a, -a}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2)/((-a + x)^2 + (y - yl)^2 + z^2)^(5/2) - 
       E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + 
         z^2] \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2)^(3/2) + (
       3 I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2 \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2)^2 - (
       I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + 
          z^2] \[Beta]) \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2) - (
       E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2 \[Beta]^2)/((-a + x)^2 + (y - yl)^2 + z^2)^(
       3/2), {yl, +b, -b}, AccuracyGoal->aG ] + 
     NIntegrate[(
       3 E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2)/((a + x)^2 + (y - yl)^2 + z^2)^(5/2) - 
       E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + 
         z^2] \[Beta])/((a + x)^2 + (y - yl)^2 + z^2)^(3/2) + (
       3 I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2 \[Beta])/((a + x)^2 + (y - yl)^2 + z^2)^2 - (
       I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + 
          z^2] \[Beta]) \[Beta])/((a + x)^2 + (y - yl)^2 + z^2) - (
       E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
          yl)^2 \[Beta]^2)/((a + x)^2 + (y - yl)^2 + z^2)^(
       3/2), {yl, -b, +b}, AccuracyGoal->aG])
If[printCalcs, Print[ ?Ey ], Nothing ]

Ez[\[Omega]_, \[Beta]_, \[Mu]r_, \[Epsilon]r_, Irf_, a_, b_, x_, y_, 
  z_] := -(I Irf/(4 Pi \[Omega]  \[Epsilon]r \[Epsilon]0)) \
(NIntegrate[(
      3 E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
         xl) z)/((x - xl)^2 + (-b + y)^2 + z^2)^(5/2) + (
      3 I E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
         xl) z \[Beta])/((x - xl)^2 + (-b + y)^2 + z^2)^2 - (
      E^(-I Sqrt[(x - xl)^2 + (-b + y)^2 + z^2] \[Beta]) (x - 
         xl) z \[Beta]^2)/((x - xl)^2 + (-b + y)^2 + z^2)^(
      3/2), {xl, -a, +a}, AccuracyGoal->aG ] + 
    NIntegrate[(
      3 E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
         xl) z)/((x - xl)^2 + (b + y)^2 + z^2)^(5/2) + (
      3 I E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
         xl) z \[Beta])/((x - xl)^2 + (b + y)^2 + z^2)^2 - (
      E^(-I Sqrt[(x - xl)^2 + (b + y)^2 + z^2] \[Beta]) (x - 
         xl) z \[Beta]^2)/((x - xl)^2 + (b + y)^2 + z^2)^(
      3/2), {xl, +a, -a}, AccuracyGoal->aG ] + 
    NIntegrate[(
      3 E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z)/((-a + x)^2 + (y - yl)^2 + z^2)^(5/2) + (
      3 I E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z \[Beta])/((-a + x)^2 + (y - yl)^2 + z^2)^2 - (
      E^(-I Sqrt[(-a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z \[Beta]^2)/((-a + x)^2 + (y - yl)^2 + z^2)^(
      3/2), {yl, +b, -b}, AccuracyGoal->aG ] + 
    NIntegrate[(
      3 E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z)/((a + x)^2 + (y - yl)^2 + z^2)^(5/2) + (
      3 I E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z \[Beta])/((a + x)^2 + (y - yl)^2 + z^2)^2 - (
      E^(-I Sqrt[(a + x)^2 + (y - yl)^2 + z^2] \[Beta]) (y - 
         yl) z \[Beta]^2)/((a + x)^2 + (y - yl)^2 + z^2)^(
      3/2), {yl, -b, +b}, AccuracyGoal->aG] )
If[printCalcs, Print[ ?Ez ], Nothing ]

Etotal[\[Omega]_, \[Beta]_, \[Mu]r_, \[Epsilon]r_, Irf_, a_, b_, x_, 
  y_, z_] := 
 Sqrt[Abs[Ex[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
      z]]^2 + Abs[
     Ey[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
      z]]^2 ]
If[printCalcs, Print[ ?Etotal ], Nothing ]

(* Poynting vector *)
S[\[Omega]_, \[Beta]_, \[Mu]r_, \[Epsilon]r_, Irf_, a_, b_, x_, y_, 
  z_] := 0.5 Norm[
   Cross[{N[
      Ex[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]], N[Ey[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x,
        y, z]], 
     N[Ez[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]}, Conjugate[{Hx[\[Beta], Irf, a, b, x, y, z], 
      Hy[\[Beta], Irf, a, b, x, y, z], 
      Hz[\[Beta], Irf, a, b, x, y, z]}]]]
If[printCalcs, Print[ ?S ], Nothing ]

(* Constants *)
\[Epsilon]0 = 8.8541878128 10^-12
\[Mu]0 = 1.25663706212 10^-6
\[Epsilon]r = 1
\[Mu]r = 1
f = 3.5 10^6
\[Omega] = 2 Pi f
Irf = 1.0
\[Lambda] = 1/(f Sqrt[\[Mu]r \[Mu]0 \[Epsilon]r \[Epsilon]0] )
\[Beta] = 2 Pi / \[Lambda]
a = 2.054/2
b = 1.65/2
Print["-- AccuracyGoal ", aG]
Print["-- Frequency    ", f, " Hz"]
Print["-- Wavelength   ", \[Lambda] , " m"]
Print["-- RF current   ", Irf , " Arms"]
Print["-- Loop width   ", a , " m"]
Print["-- Loop height  ", b , " m"]

(* Loop center *)
Print["-- Loop center magnetic field [A/m]"]
Print[ Htotal[\[Beta], Irf, a, b, 0, 0, 0] ]
Print["-- Loop center electric field [V/m]"]
Print[ Etotal[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, 0, 0, 0] ]

(* Magnetic field vector plot *)
If[ calcHfield, Print["-- Calculate magnetic vector field"], Nothing ]
If[ calcHfield, HRexyz = Table[{{x, y, z}, 
   {Re[Hx[\[Beta], Irf, a, b, x, y, z]], 
    Re[Hy[\[Beta], Irf, a, b, x, y, z]], 
    Re[Hz[\[Beta], Irf, a, b, x, y, z]]}}, {x, 0, 2, 1}, {y, 0, 2, 
   1}, {z, 0.1, 3, 1}], Nothing ]
If[ calcHfield, HImxyz = Table[{{x, y, z}, 
    {Im[Hx[\[Beta], Irf, a, b, x, y, z]], 
    Im[Hy[\[Beta], Irf, a, b, x, y, z]], 
    Im[Hz[\[Beta], Irf, a, b, x, y, z]]}}, {x, 0, 2, 1}, {y, 0, 2, 
   1}, {z, 0.1, 3, 1}], Nothing ]

If[ calcHfield, Export["/tmp/H_Real.pdf", 
  ListVectorPlot3D[HRexyz, VectorColorFunction -> "Rainbow", 
  PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]], Nothing ]
 
If[ calcHfield, Print["-- New file /tmp/H_Real.pdf"], Nothing ]

If[ calcHfield, Export["/tmp/H_Imaginary.pdf", 
  ListVectorPlot3D[HImxyz, VectorColorFunction -> "Rainbow", 
  PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]], Nothing ] 
If[ calcHfield, Print["-- New file /tmp/H_Imaginary.pdf"], Nothing ]

(* Electric field vector plot *)
Print["-- Calculate electric vector field"]
ERexyz = Table[{{x, y, z}, {Re[N[
      Ex[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], Re[
     N[Ey[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], 0}}, {x, 0, 2, 1}, {y, 0, 2, 1}, {z, 0.1, 3, 1}] 

EImxyz = Table[{{x, y, 
    z}, {Im[N[
      Ex[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], Im[
     N[Ey[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], 0}}, {x, 0, 2, 1}, {y, 0, 2, 1}, {z, 0.1, 3, 1}] 

Export["/tmp/E_Read.pdf", 
 ListVectorPlot3D[ERexyz, VectorColorFunction -> "Rainbow", 
 PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]] 
Print["-- New file /tmp/E_Real.pdf"]

EImxyz = Table[{{x, y, 
    z}, {Im[N[
      Ex[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], Im[
     N[Ey[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, y, 
       z]]], 0}}, {x, 0, 2, 1}, {y, 0, 2, 1}, {z, 0.1, 3, 1}] 

Export["/tmp/E_Imaginary.pdf",
 ListVectorPlot3D[EImxyz, VectorColorFunction -> "Rainbow", 
 PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]] 
Print["-- New file /tmp/E_Imaginary.pdf"]

(* Poynting field vector plot *)
Print["-- Calculate Poynting vector field"]
Sxyz = Table[{{1 (i - 1), 1 (j - 1), 1 (k - 1) + 0.1}, 
   N[0.5 Cross[
      ERexyz[[i]][[j]][[k]][[2]] + I EImxyz[[i]][[j]][[k]][[2]], 
      Conjugate[
       HRexyz[[i]][[j]][[k]][[2]] + 
        I HImxyz[[i]][[j]][[k]][[2]]]]]}, {i, 1, 3}, {j, 1, 3}, {k, 1,
    3}] 
SxyzIm = Table[{Sxyz[[i]][[j]][[k]][[1]], 
   Im[Sxyz[[i]][[j]][[k]][[2]]]}, {i, 1, 3}, {j, 1, 3}, {k, 1, 3}] 

Export["/tmp/S_Real.pdf", 
  ListVectorPlot3D[Re[Sxyz], VectorColorFunction -> "Rainbow", 
  PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]] 
Print["-- New file /tmp/S_Real.pdf"]

Export["/tmp/S_Imaginary.pdf",
  ListVectorPlot3D[SxyzIm, VectorColorFunction -> "Rainbow", 
  PlotLegends -> Automatic, 
 AxesLabel -> {"x[m]", "y[m]", "z[m]"}]] 
Print["-- New file /tmp/S_Imaginary.pdf"]

(* Plots *)
Print["-- Calculate magnetic field 0.1 - 6 m from loop center"]
Htotalxy0 = Table[{z, Htotal[\[Beta], Irf, a, b, x, 0, z]}, {x, 0, 2, 0.5}, {z, 0.1, 6, 0.1}]

Htotalxymb = Table[{z, Htotal[\[Beta], Irf, a, b, x, -b/2, z]}, {x, 0, 2, 0.5}, {z, 0.1, 6, 0.1}]

Export["/tmp/MagneticField.pdf", ListLinePlot[Join[Htotalxy0, Htotalxymb], 
 PlotLegends -> {"H[A/m], x = 0 m, y = 0 m", 
    "H[A/m], x = 0.5 m, y = 0 m", 
    "H[A/m], x = 1.0 m, y = 0 m", 
    "H[A/m], x = 1.5 m, y = 0 m", 
    "H[A/m], x = 2.0 m, y = 0 m", 
    "H[A/m], x = 0 m, y = -b/2 m", 
    "H[A/m], x = 0.5 m, y = -b/2 m", 
    "H[A/m], x = 1.0 m, y = -b/2 m", 
    "H[A/m], x = 1.5 m, y = -b/2 m", 
    "H[A/m], x = 2.0 m, y = -b/2 m"}, 
 AxesLabel -> {"z[m]", "A/m"}]]
Print["-- New file /tmp/MagneticField.pdf"]

Print["-- Calculate electric field 0.1 - 6 m from loop center"]
Etotalxy0 = 
 Table[{z, 
   N[Etotal[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, 0, 
     z]]}, {x, 0, 2, 0.5}, {z, 0.1, 6, 1}]

Etotalxymb = 
 Table[{z, 
   N[Etotal[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, 
     x, -b/2, z]]}, {x, 0, 2, 0.5}, {z, 0.1, 6, 1}]

Export["/tmp/ElectricField.pdf", ListLinePlot[Join[Etotalxy0, Etotalxymb], 
 PlotLegends -> {"E[V/m], x = 0 m, y = 0 m", 
   "E[V/m], x = 0.5 m, y = 0 m", 
   "E[V/m], x = 1.0 m, y = 0 m", 
   "E[V/m], x = 1.5 m, y = 0 m", 
   "E[V/m], x = 2.0 m, y = 0 m", 
   "E[V/m], x = 0 m, y = -b/2 m", 
   "E[V/m], x = 0.5 m, y = -b/2 m", 
   "E[V/m], x = 1.0 m, y = -b/2 m", 
   "E[V/m], x = 1.5 m, y = -b/2 m", 
   "E[V/m], x = 2.0 m, y = -b/2 m"}, 
 AxesLabel -> {"z[m]", "V/m"}]]
Print["-- New file /tmp/ElectricField.pdf"]

Print["-- Calculate power density 0.1 - 6 m from loop center"]
Sxy0 = Table[{z, 
   S[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, 0, z]}, {x,
    0, 2, 0.5}, {z, 0.1, 6, 1}]

Sxymb = Table[{z, 
   S[\[Omega], \[Beta], \[Mu]r, \[Epsilon]r, Irf, a, b, x, -b/2, 
    z]}, {x, 0, 2, 0.5}, {z, 0.1, 6, 1}]

Export["/tmp/PowerDensity.pdf", ListLinePlot[Join[Sxy0, Sxymb], 
 PlotLegends -> {"S[W/m^2], x = 0 m, y = 0 m", 
   "S[W/m^2], x = 0.5 m, y = 0 m", 
   "S[W/m^2], x = 1.0 m, y = 0 m", 
   "S[W/m^2], x = 1.5 m, y = 0 m", 
   "S[W/m^2], x = 2.0 m, y = 0 m", 
   "S[W/m^2], x = 0 m, y = -b/2 m", 
   "S[W/m^2], x = 0.5 m, y = -b/2 m", 
   "S[W/m^2], x = 1.0 m, y = -b/2 m", 
   "S[W/m^2], x = 1.5 m, y = -b/2 m", 
   "S[W/m^2], x = 2.0 m, y = -b/2 m"}, 
 AxesLabel -> {"z[m]", "W/m^2"}]]
Print["-- New file /tmp/PowerDensity.pdf"]
 
