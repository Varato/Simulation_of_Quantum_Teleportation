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
draw (-1.5u,-1u)--(-1.5u,7u)--(15.5u,7u)--(15.5u,-1u)--cycle withcolor 1 white;

% down line
draw (0,0)--(10u,0);
draw (11u,0)--(12u,0);
draw (13u,0)--(14u,0);
label.lft(btex $|\psi\rangle_2$ etex, (0,0));
label.rt(btex $|\phi\rangle_0$ etex, (14u,0));

% middle line
draw (0,3u)--(1u,3u);
draw (2u,3u)--(8u,3u);
drawarrow (9u,3u)--(10.5u,3u)--(10.5u,0.5u) dashed evenly;
label.top(btex $m'$ etex, (9.5u,3u));
label.lft(btex $|\psi\rangle_1$ etex, (0,3u));

% upper line
draw (0,6u)--(6u,6u);
draw (7u,6u)--(8u,6u);
drawarrow (9u,6u)--(12.5u,6u)--(12.5u, 0.5u) dashed evenly;
label.top(btex $m$ etex, (9.5u,6u));
label.lft(btex $|\phi\rangle_0$ etex, (0,6u));

% CNOT1
draw (3u,3u)--(3u,0);
draw (3u,3u) withpen pencircle scaled 8;
draw (3u,0)+(0.25u,0.25u)--(3u,0)-(0.25u,0.25u);
draw (3u,0)+(-0.25u,0.25u)--(3u,0)+(0.25u,-0.25u);

% CNOT2
draw (5u,6u)--(5u,3u);
draw (5u,6u) withpen pencircle scaled 8;
draw (5u,3u)+(0.25u,0.25u)--(5u,3u)-(0.25u,0.25u);
draw (5u,3u)+(-0.25u,0.25u)--(5u,3u)+(0.25u,-0.25u);


% Hadamard1
draw (1u,2.5u)--(1u,3.5u)--(2u,3.5u)--(2u,2.5u)--cycle withpen pencircle scaled 2;
label(btex \textbf{H} etex, (1.5u,3u));

% Hadamard2
draw (6u,5.5u)--(6u,6.5u)--(7u,6.5u)--(7u,5.5u)--cycle withpen pencircle scaled 2;
label(btex \textbf{H} etex, (6.5u,6u));

% measure1
draw (8u, 5.5u)--(8u,6.5u)--(9u,6.5u)--(9u,5.5u)--cycle withpen pencircle scaled 2;
draw halfcircle scaled 0.7u shifted (8.5u,5.9u);
drawarrow (8.5u, 5.9u)--(8.5u,5.9u)+(0.3u,0.3u);

% measure2
draw (8u, 5.5u)--(8u,6.5u)--(9u,6.5u)--(9u,5.5u)--cycle withpen pencircle scaled 2 shifted (0,-3u);
draw halfcircle scaled 0.7u shifted (8.5u,5.9u) shifted (0,-3u);
drawarrow (8.5u, 5.9u) shifted (0,-3u)--(8.5u,5.9u)+(0.3u,0.3u) shifted (0,-3u);

% X
draw (10u,-0.5u)--(10u,0.5u)--(11u,0.5u)--(11u,-0.5u)--cycle withpen pencircle scaled 2;
label(btex \textbf{X}$^{m'}$ etex, (10.5u,0));

% Z
draw (12u,-0.5u)--(12u,0.5u)--(13u,0.5u)--(13u,-0.5u)--cycle withpen pencircle scaled 2;
label(btex \textbf{Z}$^{m}$ etex, (12.5u,0));

label(btex Alice etex, (0.5u, 4.5u));
label.top(btex Bob etex, (0.5u,0u));

draw halfcircle scaled 3.4u rotated 90 shifted (0.5u,4.5u) xscaled 0.01u shifted(0.3u,0) withcolor .7 white;

endfig;

end
















