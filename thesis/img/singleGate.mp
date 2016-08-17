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

draw unitsquare scaled u withpen pencircle scaled 2;
draw (-0.8u,0.5u)--(0,0.5u);
draw (u,0.5u)--(1.8u,0.5u);
label(btex \textbf{H} etex ,(0.5u,0.5u));
label.top(btex in etex,(-0.8u,0.5u));
label.top(btex out etex,(1.8u,0.5u));

draw (-1.3u,-0.5u)--(-1.3u,1.5u)--(2.3u,1.5u)--(2.3u,-0.5u)--cycle withcolor 1 white;

endfig;

end