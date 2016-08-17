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
beginfig(1)
% white fram
draw (-7u,-1u)--(-7u,4u)--(7u,4u)--(7u,-1u)--cycle withcolor 1 white;

% upper line
draw (-3u,3u)--(1u,3u);
draw (2u,3u)--(3u,3u);

% down line
draw (-3u,0)--(3u,0);

% CNOT
draw (0,3u)--(0,0);
draw (0,3u) withpen pencircle scaled 8;
draw (0,0)+(0.25u,0.25u)--(0,0)-(0.25u,0.25u);
draw (0,0)+(-0.25u,0.25u)--(0,0)+(0.25u,-0.25u);

% Hadamard
draw (2u,2.5u)--(2u,3.5u)--(1u,3.5u)--(1u,2.5u)--cycle withpen pencircle scaled 2;
label(btex \textbf{H} etex, (1.5u,3u));

label(btex $|\beta_{xy}\rangle$etex, (-1.5u, 1.5u));
label.rt(btex $|x\rangle$ etex, (3u,3u));
label.rt(btex $|y\rangle$ etex, (3u,0));

endfig;

end