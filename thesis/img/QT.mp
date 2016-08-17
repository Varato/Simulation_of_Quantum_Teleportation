verbatimtex
%&latex
\documentclass{article}

\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{bm}

\newfont{\cyr}{wncyr10}
\begin{document}
etex


u:=25;
wi:=10;                       % width  in units u   
he:=10;                        % height in units u
hoehe:=he*u;                  % height
breite:=wi*u;                 % width
beginfig(1)

draw fullcircle scaled 2u withpen pencircle scaled 2;
draw fullcircle scaled 2u shifted (10u,0) withpen pencircle scaled 2;
draw (0,0)+(u,0)--(10u,0)-(u,0);
draw (0,0)+(0,u)--(0,4u)--(10u,4u)--(10u,0)+(0,u) dashed evenly;

label.bot(btex Classical Channel etex,(5u,4u));
label.bot(btex Quantum Channel etex, (5u,0));
label(btex Alice etex, (0,0.6u));
label(btex Bob etex, (10u,0.6u));
label(btex $|\psi \rangle_1$ etex,(0,-0.2u));
label(btex $|\psi \rangle_2$ etex,(10u,-0.2u));
label(btex $|\phi \rangle_0$ etex,(-1.5u,-1.5u));
label(btex EPR Pair etex,(5u,-3u));
draw (4u,-3u)--(0,0)+(u,-0.8u) dashed evenly;
draw (6u,-3u)--(10u,0)+(-u,-0.8u) dashed evenly;
drawarrow (-1.3u,-1.3u)--(-0.8u,-0.8u);
draw (-2u,-3.5u)--(-2u,4.5u)--(12u,4.5u)--(12u,-3.5u)--cycle withcolor 1 white;

endfig;

end