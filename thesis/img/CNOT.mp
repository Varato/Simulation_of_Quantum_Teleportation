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
draw (-3u,3u)--(3u,3u);
draw (-3u,0)--(3u,0);
draw (0,3u)--(0,0);
draw (0,3u) withpen pencircle scaled 8;
draw (0,0)+(0.25u,0.25u)--(0,0)-(0.25u,0.25u);
draw (0,0)+(-0.25u,0.25u)--(0,0)+(0.25u,-0.25u);
draw (-7u,-1u)--(-7u,4u)--(7u,4u)--(7u,-1u)--cycle withcolor 1 white;
label.lft(btex $|a\rangle=|0\rangle$ etex, (-3u,3u));
label.lft(btex $|x\rangle=a|0\rangle+b|1\rangle$ etex, (-3u,0));
label.rt(btex $|a'\rangle=|0\rangle$ etex, (3u,3u));
label.rt(btex $|x'\rangle=a|0\rangle+b|1\rangle$ etex, (3u,0));
endfig;

beginfig(2)
draw (-3u,3u)--(3u,3u);
draw (-3u,0)--(3u,0);
draw (0,3u)--(0,0);
draw (0,3u) withpen pencircle scaled 8;
draw (0,0)+(0.25u,0.25u)--(0,0)-(0.25u,0.25u);
draw (0,0)+(-0.25u,0.25u)--(0,0)+(0.25u,-0.25u);
draw (-7u,-1u)--(-7u,4u)--(7u,4u)--(7u,-1u)--cycle withcolor 1 white;
label.lft(btex $|a\rangle=|1\rangle$ etex, (-3u,3u));
label.lft(btex $|x\rangle=a|0\rangle+b|1\rangle$ etex, (-3u,0));
label.rt(btex $|a'\rangle=|1\rangle$ etex, (3u,3u));
label.rt(btex $|x'\rangle=b|0\rangle+a|1\rangle$ etex, (3u,0));
endfig;

beginfig(3)
draw (-3u,3u)--(3u,3u);
draw (-3u,0)--(3u,0);
draw (0,3u)--(0,0);
draw (0,3u) withpen pencircle scaled 8;
draw (0,0)+(0.25u,0.25u)--(0,0)-(0.25u,0.25u);
draw (0,0)+(-0.25u,0.25u)--(0,0)+(0.25u,-0.25u);
draw (-7u,-1u)--(-7u,4u)--(7u,4u)--(7u,-1u)--cycle withcolor 1 white;
label.lft(btex $\alpha|0\rangle+\beta|1\rangle$ etex, (-3u,3u));
label.lft(btex $x|0\rangle+y|1\rangle$ etex, (-3u,0));
endfig;
end