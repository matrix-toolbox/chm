function Y = A36_BH(P)
% -----------------------------------------------------------------------------
% 2023-06-12
% w.bruzda[at]uj.edu.pl
% -----------------------------------------------------------------------------
% An affine family of 2-Unitary CHM with 19 parameters stemming from BH(36, 6).
% Special (simplified) case of the general formula described in "A36.m".
% Generic defect = 61.
% -----------------------------------------------------------------------------
% Usage:
%
% >> Y = A36_BH(p)
%    where p is a vector of ANY phases in R^19
%
% >> Y = A36_BH(rand(1, 19));
%
% in particular
% >> B = A36([0 0 1/3 0 0], zeros(1, 19)); % is BH(36, 6)
% -----------------------------------------------------------------------------
% It is not known but rather unlikely that this family contains a real
% representative of CHM(36).
% -----------------------------------------------------------------------------
% Graphical representation of BH(36, 6) is shown in file: "B36_6.png".
% -----------------------------------------------------------------------------

	a = P(1);
	b = P(2);
	c = P(3);
	d = P(4);
	e = P(5);
	f = P(6);
	g = P(7);
	h = P(8);
	i = P(9);
	j = P(10);
	k = P(11);
	l = P(12);
	m = P(13);
	n = P(14);
	o = P(15);
	p = P(16);
	q = P(17);
	r = P(18);
	s = P(19);

	R = [
		a+c+f+j-l+n a+d+g+j+n a+e+n a+f+j-l+n a+g+j+n a+n   c+f+j-l+n d+g+j+n e+n f+j-l+n g+j+n +n   b+c+f+j-l+n b+d+g+j+n b+e+n b+f+j-l+n b+g+j+n b+n   a+c+f+j-l+n a+d+g+j+n a+e+n a+f+j-l+n a+g+j+n a+n   c+f+j-l+n d+g+j+n e+n f+j-l+n g+j+n +n   b+c+f+j-l+n b+d+g+j+n b+e+n b+f+j-l+n b+g+j+n b+n;
		a+c+f+k-m+o a+d+g+k+o a+e+o a+f+k-m+o a+g+k+o a+o   c+f+k-m+o d+g+k+o e+o f+k-m+o g+k+o +o   b+c+f+k-m+o b+d+g+k+o b+e+o b+f+k-m+o b+g+k+o b+o   a+c+f+k-m+o a+d+g+k+o a+e+o a+f+k-m+o a+g+k+o a+o   c+f+k-m+o d+g+k+o e+o f+k-m+o g+k+o +o   b+c+f+k-m+o b+d+g+k+o b+e+o b+f+k-m+o b+g+k+o b+o;
		a+c+f+j-l+p a+d+g+j+p a+e+p a+f+j-l+p a+g+j+p a+p   c+f+j-l+p d+g+j+p e+p f+j-l+p g+j+p +p   b+c+f+j-l+p b+d+g+j+p b+e+p b+f+j-l+p b+g+j+p b+p   a+c+f+j-l+p a+d+g+j+p a+e+p a+f+j-l+p a+g+j+p a+p   c+f+j-l+p d+g+j+p e+p f+j-l+p g+j+p +p   b+c+f+j-l+p b+d+g+j+p b+e+p b+f+j-l+p b+g+j+p b+p;
		a+c+f+k-m+q a+d+g+k+q a+e+q a+f+k-m+q a+g+k+q a+q   c+f+k-m+q d+g+k+q e+q f+k-m+q g+k+q +q   b+c+f+k-m+q b+d+g+k+q b+e+q b+f+k-m+q b+g+k+q b+q   a+c+f+k-m+q a+d+g+k+q a+e+q a+f+k-m+q a+g+k+q a+q   c+f+k-m+q d+g+k+q e+q f+k-m+q g+k+q +q   b+c+f+k-m+q b+d+g+k+q b+e+q b+f+k-m+q b+g+k+q b+q;
		a+c+f+j-l+r a+d+g+j+r a+e+r a+f+j-l+r a+g+j+r a+r   c+f+j-l+r d+g+j+r e+r f+j-l+r g+j+r +r   b+c+f+j-l+r b+d+g+j+r b+e+r b+f+j-l+r b+g+j+r b+r   a+c+f+j-l+r a+d+g+j+r a+e+r a+f+j-l+r a+g+j+r a+r   c+f+j-l+r d+g+j+r e+r f+j-l+r g+j+r +r   b+c+f+j-l+r b+d+g+j+r b+e+r b+f+j-l+r b+g+j+r b+r;
		a+c+f+k-m+s a+d+g+k+s a+e+s a+f+k-m+s a+g+k+s a+s   c+f+k-m+s d+g+k+s e+s f+k-m+s g+k+s +s   b+c+f+k-m+s b+d+g+k+s b+e+s b+f+k-m+s b+g+k+s b+s   a+c+f+k-m+s a+d+g+k+s a+e+s a+f+k-m+s a+g+k+s a+s   c+f+k-m+s d+g+k+s e+s f+k-m+s g+k+s +s   b+c+f+k-m+s b+d+g+k+s b+e+s b+f+k-m+s b+g+k+s b+s;

		a+c+j+n a+d+n a+e+j-l+n a+j+n a+n a+j-l+n   c+j+n d+n e+j-l+n j+n, n, j-l+n   b+c+j+n b+d+n b+e+j-l+n b+j+n b+n b+j-l+n   a+c+j+n a+d+n a+e+j-l+n a+j+n a+n a+j-l+n   c+j+n d+n e+j-l+n +j+n +n +j-l+n   b+c+j+n b+d+n b+e+j-l+n b+j+n b+n b+j-l+n;
		a+c+k+o a+d+o a+e+k-m+o a+k+o a+o a+k-m+o   c+k+o d+o e+k-m+o k+o, o, k-m+o   b+c+k+o b+d+o b+e+k-m+o b+k+o b+o b+k-m+o   a+c+k+o a+d+o a+e+k-m+o a+k+o a+o a+k-m+o   c+k+o d+o e+k-m+o +k+o +o +k-m+o   b+c+k+o b+d+o b+e+k-m+o b+k+o b+o b+k-m+o;
		a+c+j+p a+d+p a+e+j-l+p a+j+p a+p a+j-l+p   c+j+p d+p e+j-l+p j+p, p, j-l+p   b+c+j+p b+d+p b+e+j-l+p b+j+p b+p b+j-l+p   a+c+j+p a+d+p a+e+j-l+p a+j+p a+p a+j-l+p   c+j+p d+p e+j-l+p +j+p +p +j-l+p   b+c+j+p b+d+p b+e+j-l+p b+j+p b+p b+j-l+p;
		a+c+k+q a+d+q a+e+k-m+q a+k+q a+q a+k-m+q   c+k+q d+q e+k-m+q k+q, q, k-m+q   b+c+k+q b+d+q b+e+k-m+q b+k+q b+q b+k-m+q   a+c+k+q a+d+q a+e+k-m+q a+k+q a+q a+k-m+q   c+k+q d+q e+k-m+q +k+q +q +k-m+q   b+c+k+q b+d+q b+e+k-m+q b+k+q b+q b+k-m+q;
		a+c+j+r a+d+r a+e+j-l+r a+j+r a+r a+j-l+r   c+j+r d+r e+j-l+r j+r, r, j-l+r   b+c+j+r b+d+r b+e+j-l+r b+j+r b+r b+j-l+r   a+c+j+r a+d+r a+e+j-l+r a+j+r a+r a+j-l+r   c+j+r d+r e+j-l+r +j+r +r +j-l+r   b+c+j+r b+d+r b+e+j-l+r b+j+r b+r b+j-l+r;
		a+c+k+s a+d+s a+e+k-m+s a+k+s a+s a+k-m+s   c+k+s d+s e+k-m+s k+s, s, k-m+s   b+c+k+s b+d+s b+e+k-m+s b+k+s b+s b+k-m+s   a+c+k+s a+d+s a+e+k-m+s a+k+s a+s a+k-m+s   c+k+s d+s e+k-m+s +k+s +s +k-m+s   b+c+k+s b+d+s b+e+k-m+s b+k+s b+s b+k-m+s;

		a+c+i+n a+d+h+j-l+n a+e+j+n a+i+n a+h+j-l+n a+j+n   c+i+n d+h+j-l+n e+j+n i+n h+j-l+n +j+n   b+c+i+n b+d+h+j-l+n b+e+j+n b+i+n b+h+j-l+n b+j+n   a+c+i+n a+d+h+j-l+n a+e+j+n a+i+n a+h+j-l+n a+j+n   c+i+n d+h+j-l+n e+j+n i+n h+j-l+n +j+n   b+c+i+n b+d+h+j-l+n b+e+j+n b+i+n b+h+j-l+n b+j+n;
		a+c+i+o a+d+h+k-m+o a+e+k+o a+i+o a+h+k-m+o a+k+o   c+i+o d+h+k-m+o e+k+o i+o h+k-m+o +k+o   b+c+i+o b+d+h+k-m+o b+e+k+o b+i+o b+h+k-m+o b+k+o   a+c+i+o a+d+h+k-m+o a+e+k+o a+i+o a+h+k-m+o a+k+o   c+i+o d+h+k-m+o e+k+o i+o h+k-m+o +k+o   b+c+i+o b+d+h+k-m+o b+e+k+o b+i+o b+h+k-m+o b+k+o;
		a+c+i+p a+d+h+j-l+p a+e+j+p a+i+p a+h+j-l+p a+j+p   c+i+p d+h+j-l+p e+j+p i+p h+j-l+p +j+p   b+c+i+p b+d+h+j-l+p b+e+j+p b+i+p b+h+j-l+p b+j+p   a+c+i+p a+d+h+j-l+p a+e+j+p a+i+p a+h+j-l+p a+j+p   c+i+p d+h+j-l+p e+j+p i+p h+j-l+p +j+p   b+c+i+p b+d+h+j-l+p b+e+j+p b+i+p b+h+j-l+p b+j+p;
		a+c+i+q a+d+h+k-m+q a+e+k+q a+i+q a+h+k-m+q a+k+q   c+i+q d+h+k-m+q e+k+q i+q h+k-m+q +k+q   b+c+i+q b+d+h+k-m+q b+e+k+q b+i+q b+h+k-m+q b+k+q   a+c+i+q a+d+h+k-m+q a+e+k+q a+i+q a+h+k-m+q a+k+q   c+i+q d+h+k-m+q e+k+q i+q h+k-m+q +k+q   b+c+i+q b+d+h+k-m+q b+e+k+q b+i+q b+h+k-m+q b+k+q;
		a+c+i+r a+d+h+j-l+r a+e+j+r a+i+r a+h+j-l+r a+j+r   c+i+r d+h+j-l+r e+j+r i+r h+j-l+r +j+r   b+c+i+r b+d+h+j-l+r b+e+j+r b+i+r b+h+j-l+r b+j+r   a+c+i+r a+d+h+j-l+r a+e+j+r a+i+r a+h+j-l+r a+j+r   c+i+r d+h+j-l+r e+j+r i+r h+j-l+r +j+r   b+c+i+r b+d+h+j-l+r b+e+j+r b+i+r b+h+j-l+r b+j+r;
		a+c+i+s a+d+h+k-m+s a+e+k+s a+i+s a+h+k-m+s a+k+s   c+i+s d+h+k-m+s e+k+s i+s h+k-m+s +k+s   b+c+i+s b+d+h+k-m+s b+e+k+s b+i+s b+h+k-m+s b+k+s   a+c+i+s a+d+h+k-m+s a+e+k+s a+i+s a+h+k-m+s a+k+s   c+i+s d+h+k-m+s e+k+s i+s h+k-m+s +k+s   b+c+i+s b+d+h+k-m+s b+e+k+s b+i+s b+h+k-m+s b+k+s;

		a+c+f+j-l+n a+d+g+j+n a+e+n a+f+j-l+n a+g+j+n a+n   c+f+j-l+n d+g+j+n e+n f+j-l+n g+j+n +n   b+c+f+j-l+n b+d+g+j+n b+e+n b+f+j-l+n b+g+j+n b+n   a+c+f+j-l+n a+d+g+j+n a+e+n a+f+j-l+n a+g+j+n a+n   c+f+j-l+n d+g+j+n e+n f+j-l+n g+j+n +n   b+c+f+j-l+n b+d+g+j+n b+e+n b+f+j-l+n b+g+j+n b+n;
		a+c+f+k-m+o a+d+g+k+o a+e+o a+f+k-m+o a+g+k+o a+o   c+f+k-m+o d+g+k+o e+o f+k-m+o g+k+o +o   b+c+f+k-m+o b+d+g+k+o b+e+o b+f+k-m+o b+g+k+o b+o   a+c+f+k-m+o a+d+g+k+o a+e+o a+f+k-m+o a+g+k+o a+o   c+f+k-m+o d+g+k+o e+o f+k-m+o g+k+o +o   b+c+f+k-m+o b+d+g+k+o b+e+o b+f+k-m+o b+g+k+o b+o;
		a+c+f+j-l+p a+d+g+j+p a+e+p a+f+j-l+p a+g+j+p a+p   c+f+j-l+p d+g+j+p e+p f+j-l+p g+j+p +p   b+c+f+j-l+p b+d+g+j+p b+e+p b+f+j-l+p b+g+j+p b+p   a+c+f+j-l+p a+d+g+j+p a+e+p a+f+j-l+p a+g+j+p a+p   c+f+j-l+p d+g+j+p e+p f+j-l+p g+j+p +p   b+c+f+j-l+p b+d+g+j+p b+e+p b+f+j-l+p b+g+j+p b+p;
		a+c+f+k-m+q a+d+g+k+q a+e+q a+f+k-m+q a+g+k+q a+q   c+f+k-m+q d+g+k+q e+q f+k-m+q g+k+q +q   b+c+f+k-m+q b+d+g+k+q b+e+q b+f+k-m+q b+g+k+q b+q   a+c+f+k-m+q a+d+g+k+q a+e+q a+f+k-m+q a+g+k+q a+q   c+f+k-m+q d+g+k+q e+q f+k-m+q g+k+q +q   b+c+f+k-m+q b+d+g+k+q b+e+q b+f+k-m+q b+g+k+q b+q;
		a+c+f+j-l+r a+d+g+j+r a+e+r a+f+j-l+r a+g+j+r a+r   c+f+j-l+r d+g+j+r e+r f+j-l+r g+j+r +r   b+c+f+j-l+r b+d+g+j+r b+e+r b+f+j-l+r b+g+j+r b+r   a+c+f+j-l+r a+d+g+j+r a+e+r a+f+j-l+r a+g+j+r a+r   c+f+j-l+r d+g+j+r e+r f+j-l+r g+j+r +r   b+c+f+j-l+r b+d+g+j+r b+e+r b+f+j-l+r b+g+j+r b+r;
		a+c+f+k-m+s a+d+g+k+s a+e+s a+f+k-m+s a+g+k+s a+s   c+f+k-m+s d+g+k+s e+s f+k-m+s g+k+s +s   b+c+f+k-m+s b+d+g+k+s b+e+s b+f+k-m+s b+g+k+s b+s   a+c+f+k-m+s a+d+g+k+s a+e+s a+f+k-m+s a+g+k+s a+s   c+f+k-m+s d+g+k+s e+s f+k-m+s g+k+s +s   b+c+f+k-m+s b+d+g+k+s b+e+s b+f+k-m+s b+g+k+s b+s;

		a+c+j+n a+d+n a+e+j-l+n a+j+n a+n a+j-l+n   c+j+n d+n e+j-l+n +j+n +n +j-l+n   b+c+j+n b+d+n b+e+j-l+n b+j+n b+n b+j-l+n   a+c+j+n a+d+n a+e+j-l+n a+j+n a+n a+j-l+n   c+j+n d+n e+j-l+n +j+n +n +j-l+n   b+c+j+n b+d+n b+e+j-l+n b+j+n b+n b+j-l+n;
		a+c+k+o a+d+o a+e+k-m+o a+k+o a+o a+k-m+o   c+k+o d+o e+k-m+o +k+o +o +k-m+o   b+c+k+o b+d+o b+e+k-m+o b+k+o b+o b+k-m+o   a+c+k+o a+d+o a+e+k-m+o a+k+o a+o a+k-m+o   c+k+o d+o e+k-m+o +k+o +o +k-m+o   b+c+k+o b+d+o b+e+k-m+o b+k+o b+o b+k-m+o;
		a+c+j+p a+d+p a+e+j-l+p a+j+p a+p a+j-l+p   c+j+p d+p e+j-l+p +j+p +p +j-l+p   b+c+j+p b+d+p b+e+j-l+p b+j+p b+p b+j-l+p   a+c+j+p a+d+p a+e+j-l+p a+j+p a+p a+j-l+p   c+j+p d+p e+j-l+p +j+p +p +j-l+p   b+c+j+p b+d+p b+e+j-l+p b+j+p b+p b+j-l+p;
		a+c+k+q a+d+q a+e+k-m+q a+k+q a+q a+k-m+q   c+k+q d+q e+k-m+q +k+q +q +k-m+q   b+c+k+q b+d+q b+e+k-m+q b+k+q b+q b+k-m+q   a+c+k+q a+d+q a+e+k-m+q a+k+q a+q a+k-m+q   c+k+q d+q e+k-m+q +k+q +q +k-m+q   b+c+k+q b+d+q b+e+k-m+q b+k+q b+q b+k-m+q;
		a+c+j+r a+d+r a+e+j-l+r a+j+r a+r a+j-l+r   c+j+r d+r e+j-l+r +j+r +r +j-l+r   b+c+j+r b+d+r b+e+j-l+r b+j+r b+r b+j-l+r   a+c+j+r a+d+r a+e+j-l+r a+j+r a+r a+j-l+r   c+j+r d+r e+j-l+r +j+r +r +j-l+r   b+c+j+r b+d+r b+e+j-l+r b+j+r b+r b+j-l+r;
		a+c+k+s a+d+s a+e+k-m+s a+k+s a+s a+k-m+s   c+k+s d+s e+k-m+s +k+s +s +k-m+s   b+c+k+s b+d+s b+e+k-m+s b+k+s b+s b+k-m+s   a+c+k+s a+d+s a+e+k-m+s a+k+s a+s a+k-m+s   c+k+s d+s e+k-m+s +k+s +s +k-m+s   b+c+k+s b+d+s b+e+k-m+s b+k+s b+s b+k-m+s;

		a+c+i+n a+d+h+j-l+n a+e+j+n a+i+n a+h+j-l+n a+j+n   c+i+n d+h+j-l+n e+j+n i+n h+j-l+n +j+n   b+c+i+n b+d+h+j-l+n b+e+j+n b+i+n b+h+j-l+n b+j+n   a+c+i+n a+d+h+j-l+n a+e+j+n a+i+n a+h+j-l+n a+j+n   c+i+n d+h+j-l+n e+j+n i+n h+j-l+n +j+n   b+c+i+n b+d+h+j-l+n b+e+j+n b+i+n b+h+j-l+n b+j+n;
		a+c+i+o a+d+h+k-m+o a+e+k+o a+i+o a+h+k-m+o a+k+o   c+i+o d+h+k-m+o e+k+o i+o h+k-m+o +k+o   b+c+i+o b+d+h+k-m+o b+e+k+o b+i+o b+h+k-m+o b+k+o   a+c+i+o a+d+h+k-m+o a+e+k+o a+i+o a+h+k-m+o a+k+o   c+i+o d+h+k-m+o e+k+o i+o h+k-m+o +k+o   b+c+i+o b+d+h+k-m+o b+e+k+o b+i+o b+h+k-m+o b+k+o;
		a+c+i+p a+d+h+j-l+p a+e+j+p a+i+p a+h+j-l+p a+j+p   c+i+p d+h+j-l+p e+j+p i+p h+j-l+p +j+p   b+c+i+p b+d+h+j-l+p b+e+j+p b+i+p b+h+j-l+p b+j+p   a+c+i+p a+d+h+j-l+p a+e+j+p a+i+p a+h+j-l+p a+j+p   c+i+p d+h+j-l+p e+j+p i+p h+j-l+p +j+p   b+c+i+p b+d+h+j-l+p b+e+j+p b+i+p b+h+j-l+p b+j+p;
		a+c+i+q a+d+h+k-m+q a+e+k+q a+i+q a+h+k-m+q a+k+q   c+i+q d+h+k-m+q e+k+q i+q h+k-m+q +k+q   b+c+i+q b+d+h+k-m+q b+e+k+q b+i+q b+h+k-m+q b+k+q   a+c+i+q a+d+h+k-m+q a+e+k+q a+i+q a+h+k-m+q a+k+q   c+i+q d+h+k-m+q e+k+q i+q h+k-m+q +k+q   b+c+i+q b+d+h+k-m+q b+e+k+q b+i+q b+h+k-m+q b+k+q;
		a+c+i+r a+d+h+j-l+r a+e+j+r a+i+r a+h+j-l+r a+j+r   c+i+r d+h+j-l+r e+j+r i+r h+j-l+r +j+r   b+c+i+r b+d+h+j-l+r b+e+j+r b+i+r b+h+j-l+r b+j+r   a+c+i+r a+d+h+j-l+r a+e+j+r a+i+r a+h+j-l+r a+j+r   c+i+r d+h+j-l+r e+j+r i+r h+j-l+r +j+r   b+c+i+r b+d+h+j-l+r b+e+j+r b+i+r b+h+j-l+r b+j+r;
		a+c+i+s a+d+h+k-m+s a+e+k+s a+i+s a+h+k-m+s a+k+s   c+i+s d+h+k-m+s e+k+s i+s h+k-m+s +k+s   b+c+i+s b+d+h+k-m+s b+e+k+s b+i+s b+h+k-m+s b+k+s   a+c+i+s a+d+h+k-m+s a+e+k+s a+i+s a+h+k-m+s a+k+s   c+i+s d+h+k-m+s e+k+s i+s h+k-m+s +k+s   b+c+i+s b+d+h+k-m+s b+e+k+s b+i+s b+h+k-m+s b+k+s;
	];

	LA = exp(2j*pi*[
        1  5  5  5  1  3  4  2  4  4  2  4  4  2  0  0  0  2  1  5  5  5  1  3  0  4  0  0  4  0  0  0  2  4  2  0 ;
        5  5  5  1  3  1  1  3  3  1  3  3  0  4  4  4  0  2  2  2  2  4  0  4  0  2  2  0  2  2  1  3  5  3  1  1 ;
        3  3  5  1  5  3  0  0  4  0  0  4  0  0  0  2  4  2  3  3  5  1  5  3  2  2  0  2  2  0  2  4  2  0  0  0 ;
        1  3  5  3  1  1  3  1  3  3  1  3  2  2  4  0  4  2  4  0  2  0  4  4  2  0  2  2  0  2  3  1  5  5  5  1 ;
        5  1  5  3  3  3  2  4  4  2  4  4  2  4  0  4  2  2  5  1  5  3  3  3  4  0  0  4  0  0  4  2  2  2  4  0 ;
        3  1  5  5  5  1  5  5  3  5  5  3  4  0  4  2  2  2  0  4  2  2  2  4  4  4  2  4  4  2  5  5  5  1  3  1 ;
        4  2  4  1  5  1  5  3  1  4  4  0  3  1  1  4  0  2  0  4  0  3  1  3  1  1  3  2  0  4  3  1  1  4  0  2 ;
        2  4  4  5  1  1  2  0  0  3  5  1  5  5  5  4  0  4  1  3  3  4  0  0  3  5  1  2  0  0  2  2  2  1  3  1 ;
        2  2  0  5  5  3  3  3  3  2  4  2  1  1  3  2  0  4  4  4  2  1  1  5  5  1  5  0  0  0  1  1  3  2  0  4 ;
        0  4  0  3  1  3  0  0  2  1  5  3  3  5  1  2  0  0  5  3  5  2  0  2  1  5  3  0  0  2  0  2  4  5  3  3 ;
        0  2  2  3  5  5  1  3  5  0  4  4  5  1  5  0  0  0  2  4  4  5  1  1  3  1  1  4  0  2  5  1  5  0  0  0 ;
        4  4  2  1  1  5  4  0  4  5  5  5  1  5  3  0  0  2  3  3  1  0  0  4  5  5  5  4  0  4  4  2  0  3  3  5 ;
        3  3  3  5  1  5  2  2  4  0  4  2  0  0  4  0  0  4  5  1  5  3  3  3  2  2  4  0  4  2  4  4  2  4  4  2 ;
        1  1  3  5  3  1  5  1  3  1  5  5  2  0  2  2  0  2  2  0  4  4  4  0  2  4  0  4  2  2  3  1  3  3  1  3 ;
        3  5  1  5  3  3  2  4  2  0  0  0  0  2  2  0  2  2  5  3  3  3  5  1  2  4  2  0  0  0  4  0  0  4  0  0 ;
        1  3  1  5  5  5  5  3  1  1  1  3  2  2  0  2  2  0  2  2  2  4  0  4  2  0  4  4  4  0  3  3  1  3  3  1 ;
        3  1  5  5  5  1  2  0  0  0  2  4  0  4  0  0  4  0  5  5  1  3  1  5  2  0  0  0  2  4  4  2  4  4  2  4 ;
        1  5  5  5  1  3  5  5  5  1  3  1  2  4  4  2  4  4  2  4  0  4  2  2  2  2  2  4  0  4  3  5  5  3  5  5 ;
        1  2  5  2  1  0  0  1  0  3  4  3  0  3  2  1  2  3  1  2  5  2  1  0  4  5  4  1  2  1  4  5  0  3  0  5 ;
        5  2  5  4  3  4  3  2  5  0  5  2  4  3  2  3  4  1  2  5  2  1  0  1  4  3  0  1  0  3  3  4  1  4  3  2 ;
        3  0  5  4  5  0  2  5  0  5  2  3  2  1  2  3  0  3  3  0  5  4  5  0  0  3  4  3  0  1  0  3  0  5  4  5 ;
        1  0  5  0  1  4  5  0  5  2  3  2  0  1  2  5  2  1  4  3  2  3  4  1  0  1  0  3  4  3  5  2  1  0  1  2 ;
        5  4  5  0  3  0  4  3  0  1  0  3  4  5  2  5  4  3  5  4  5  0  3  0  2  1  4  5  4  1  2  1  0  1  2  5 ;
        3  4  5  2  5  4  1  4  5  4  1  2  2  5  2  1  0  1  0  1  2  5  2  1  2  5  0  5  2  3  1  0  1  2  5  2 ;
        0  1  0  0  1  0  1  4  3  5  0  1  3  4  1  1  0  5  4  5  4  4  5  4  5  0  1  1  4  3  3  4  1  1  0  5 ;
        4  3  0  4  3  0  0  5  4  2  3  0  5  2  5  1  0  1  5  4  1  5  4  1  5  0  3  3  2  1  2  5  2  4  3  4 ;
        4  1  2  4  1  2  5  4  5  3  0  3  1  4  3  5  0  1  2  5  0  2  5  0  3  0  3  5  4  5  1  4  3  5  0  1 ;
        2  3  2  2  3  2  4  5  0  0  3  2  3  2  1  5  0  3  3  4  3  3  4  3  3  0  5  1  2  3  0  5  4  2  3  0 ;
        2  1  4  2  1  4  3  4  1  1  0  5  5  4  5  3  0  3  0  5  2  0  5  2  1  0  5  3  4  1  5  4  5  3  0  3 ;
        0  3  4  0  3  4  2  5  2  4  3  4  1  2  3  3  0  5  1  4  5  1  4  5  1  0  1  5  2  5  4  5  0  0  3  2 ;
        5  4  5  0  3  0  2  5  4  3  4  5  4  1  2  1  4  5  3  0  3  2  1  2  2  5  4  3  4  5  0  3  4  3  0  1 ;
        5  0  1  4  1  0  5  4  3  4  5  2  0  1  0  3  4  3  4  1  0  5  0  1  2  1  0  1  2  5  5  0  5  2  3  2 ;
        5  0  3  0  5  4  2  1  2  3  0  3  4  3  0  1  0  3  3  2  1  2  3  0  2  1  2  3  0  3  0  5  2  3  2  5 ;
        5  2  5  4  3  4  5  0  1  4  1  0  0  3  4  3  0  1  4  3  4  5  2  5  2  3  4  1  4  3  5  2  3  2  5  0 ;
        5  2  1  0  1  2  2  3  0  3  2  1  4  5  4  1  2  1  3  4  5  2  5  4  2  3  0  3  2  1  0  1  0  3  4  3 ;
        5  4  3  4  5  2  5  2  5  4  3  4  0  5  2  3  2  5  4  5  2  5  4  3  2  5  2  1  0  1  5  4  1  2  1  4 ;
	 ] / 6);

	Y = LA .* exp(2j*pi*R/6);

end
