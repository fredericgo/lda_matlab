function v = psi(x)
	x += 6;
	v = 1 ./ (x .* x);
	v=(((0.004166666666667*v-0.003968253986254).* v +
	0.008333333333333) .* v-0.083333333333333) .* v;
	v = v + log(x)- (0.5./x) - (1 ./ (x-1)) - (1./(x-2)) - (1 ./ (x-3)) - (1 ./ (x-4)) - (1 ./ (x-5)) - (1 ./(x-6));
end
