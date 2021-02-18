% Generated using GBSolver generator v1.1 Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [w d c b meta M] = p4pf_quat101(u11, X11)

    meta = {};

	%monomials appearing in matrix M:
	%    d^3*w^2 c*d^2*w^2 c^2*d*w^2 c^3*w^2 b*d^2*w^2 b*d^3*w b*c*d*w^2 b*c*d^2*w b*c^2*w^2 b*c^2*d*w b^2*d*w^2 
	%    b^2*d^2*w b^2*c*w^2 b^2*c*d*w b^3*w^2 b^3*d*w d^2*w^2 d^3*w d^4 c*d*w^2 c*d^2*w c*d^3 
	%    c^2*w^2 c^2*d*w c^2*d^2 c^3*w c^3*d c^4 b*d*w^2 b*d^2*w b*d^3 b*c*w^2 b*c*d*w 
	%    b*c*d^2 b*c^2*w b*c^2*d b*c^3 b^2*w^2 b^2*d*w b^2*d^2 b^2*c*w b^2*c*d b^2*c^2 b^3*w 
	%    b^3*d b^3*c b^4 d*w^2 d^2*w d^3 c*w^2 c*d*w c*d^2 c^2*w c^2*d 
	%    c^3 b*w^2 b*d*w b*d^2 b*c*w b*c*d b*c^2 b^2*w b^2*d b^2*c b^3 
	%    w^2 d*w d^2 c*w c*d c^2 b*w b*d b*c b^2 w 
	%    d c b 1 
	M = zeros(65, 81);
	M([310, 374, 438, 502, 631, 695, 759, 888, 1017, 1837, 1901, 2095, 2484, 3708]) = -2;
	M([1155, 1219, 1348, 1412, 1606, 1930, 1994, 2188, 2577, 3137, 3201, 3395, 3784, 4423]) = -1;
	M([1285, 1349, 1478, 1542, 1671, 2060, 2124, 2253, 2642, 3267, 3331, 3460, 3849, 4488]) = 2;
	M([1545, 1609, 1673, 1737, 1801, 2255, 2319, 2383, 2772, 3462, 3526, 3590, 3979, 4618]) = -1;
	M([2130, 2194, 2258, 2322, 2386, 2645, 2709, 2773, 2967, 3852, 3916, 3980, 4174, 4813]) = 2/u11(2)-2*u11(1)/u11(2);
	M([2520, 2584, 2648, 2712, 2776, 2840, 2904, 2968, 3032, 4047, 4111, 4175, 4239, 4878]) = 1;
	M([3170, 3234, 3363, 3427, 3556, 3750, 3814, 3943, 4137, 4372, 4436, 4565, 4759, 5008]) = 2/u11(2)-2*u11(1)/u11(2);
	M([4405, 4469, 4533, 4597, 4661, 4725, 4789, 4853, 4917, 4957, 5021, 5085, 5149, 5203]) = 1;
	M([121, 185, 444, 508, 637, 896, 1256, 1320, 1514, 2098, 3319]) = 2*u11(6)*X11(8);
	M([316, 445, 704, 768, 897, 1026, 1841, 1905, 2099, 2488, 3709]) = 2*u11(6)*X11(7);
	M([1161, 1355, 1939, 2003, 2197, 2586, 3141, 3205, 3399, 3788, 4424]) = X11(8);
	M([1291, 1485, 2069, 2133, 2262, 2651, 3271, 3335, 3464, 3853, 4489]) = -2*u11(6)*X11(7);
	M([1551, 1680, 2264, 2328, 2392, 2781, 3466, 3530, 3594, 3983, 4619]) = -X11(8);
	M([1876, 2070, 2459, 2523, 2652, 2846, 3661, 3725, 3854, 4048, 4684]) = 2*u11(6)*X11(8);
	M([2136, 2265, 2654, 2718, 2782, 2976, 3856, 3920, 3984, 4178, 4814]) = -2*u11(6)/u11(2)+2-2*X11(7);
	M([2526, 2655, 2849, 2913, 2977, 3041, 4051, 4115, 4179, 4243, 4879]) = X11(8);
	M([3176, 3370, 3759, 3823, 3952, 4146, 4376, 4440, 4569, 4763, 5009]) = -2*u11(6)/u11(2)+2-2*X11(7);
	M([4411, 4540, 4734, 4798, 4862, 4926, 4961, 5025, 5089, 5153, 5204]) = -X11(8);
	M([59, 123, 317, 1071, 1135, 1199, 1328, 1392, 1586, 1910, 1974, 2168, 2557, 3129, 3193, 3387, 3776, 4421]) = u11(6)*X11(7)-u11(5)*X11(8);
	M([189, 253, 577, 1461, 1525, 1589, 1653, 1717, 1781, 2235, 2299, 2363, 2752, 3454, 3518, 3582, 3971, 4616]) = u11(6)*X11(7)+u11(5)*X11(8);
	M([449, 578, 837, 2046, 2110, 2174, 2238, 2302, 2366, 2625, 2689, 2753, 2947, 3844, 3908, 3972, 4166, 4811]) = 2*u11(5)*X11(7)+2*u11(6)*u11(1)/u11(2)-2*u11(6)*X11(8)-2*u11(5);
	M([709, 838, 967, 2436, 2500, 2564, 2628, 2692, 2756, 2820, 2884, 2948, 3012, 4039, 4103, 4167, 4231, 4876]) = -u11(5)*X11(8)-u11(6)*X11(7);
	M([1099, 1293, 1877, 3086, 3150, 3214, 3343, 3407, 3536, 3730, 3794, 3923, 4117, 4364, 4428, 4557, 4751, 5006]) = 2*u11(6)*X11(8)+2*u11(6)*u11(1)/u11(2)-2*u11(5)+2*u11(5)*X11(7);
	M([3114, 3308, 3697, 4321, 4385, 4449, 4513, 4577, 4641, 4705, 4769, 4833, 4897, 4949, 5013, 5077, 5141, 5201]) = -u11(6)*X11(7)+u11(5)*X11(8);
	M([63, 127, 321, 385, 1075, 1139, 1333, 1917, 3125]) = u11(8)*X11(12);
	M([128, 192, 451, 515, 1270, 1334, 1528, 2112, 3320]) = 2*u11(8)*X11(11);
	M([193, 257, 581, 645, 1465, 1529, 1658, 2242, 3450]) = -u11(8)*X11(12);
	M([323, 452, 711, 775, 1855, 1919, 2113, 2502, 3710]) = 2*u11(8)*X11(10);
	M([713, 842, 971, 1035, 2440, 2504, 2633, 2827, 4035]) = -u11(8)*X11(12);
	M([1168, 1362, 1946, 2010, 3155, 3219, 3413, 3802, 4425]) = X11(11);
	M([1298, 1492, 2076, 2140, 3285, 3349, 3478, 3867, 4490]) = -2*u11(8)*X11(10);
	M([1363, 1557, 2141, 2205, 3350, 3414, 3543, 3932, 4555]) = -2*X11(12);
	M([1558, 1687, 2271, 2335, 3480, 3544, 3608, 3997, 4620]) = -X11(11);
	M([1883, 2077, 2466, 2530, 3675, 3739, 3868, 4062, 4685]) = 2*u11(8)*X11(11);
	M([2143, 2272, 2661, 2725, 3870, 3934, 3998, 4192, 4815]) = 2-2*u11(8)/u11(2)-2*X11(10);
	M([2533, 2662, 2856, 2920, 4065, 4129, 4193, 4257, 4880]) = X11(11);
	M([3118, 3312, 3701, 3765, 4325, 4389, 4518, 4712, 4945]) = u11(8)*X11(12);
	M([3183, 3377, 3766, 3830, 4390, 4454, 4583, 4777, 5010]) = 2-2*u11(8)/u11(2)-2*X11(10);
	M([3768, 3897, 4091, 4155, 4715, 4779, 4843, 4907, 5140]) = 2*X11(12);
	M([4418, 4547, 4741, 4805, 4975, 5039, 5103, 5167, 5205]) = -X11(11);
	M([130, 324, 1081, 1405, 1599, 1988, 2182, 2571, 3133, 3197, 3391, 3780, 4422]) = u11(8)*X11(10)-u11(7)*X11(11);
	M([195, 454, 1276, 1600, 1729, 2183, 2312, 2701, 3328, 3392, 3521, 3910, 4552]) = 2*u11(7)*X11(12);
	M([260, 584, 1471, 1730, 1794, 2313, 2377, 2766, 3458, 3522, 3586, 3975, 4617]) = u11(7)*X11(11)+u11(8)*X11(10);
	M([455, 714, 1861, 2185, 2314, 2573, 2702, 2896, 3718, 3782, 3911, 4105, 4747]) = -2*u11(8)*X11(12);
	M([585, 844, 2056, 2315, 2379, 2703, 2767, 2961, 3848, 3912, 3976, 4170, 4812]) = -2*u11(8)*X11(11)+2*u11(8)*u11(1)/u11(2)-2*u11(7)+2*u11(7)*X11(10);
	M([845, 974, 2446, 2705, 2769, 2898, 2962, 3026, 4043, 4107, 4171, 4235, 4877]) = -u11(7)*X11(11)-u11(8)*X11(10);
	M([1300, 1884, 3096, 3420, 3549, 3808, 3937, 4131, 4368, 4432, 4561, 4755, 5007]) = 2*u11(7)*X11(10)-2*u11(7)+2*u11(8)*u11(1)/u11(2)+2*u11(8)*X11(11);
	M([1495, 2079, 3291, 3550, 3614, 3938, 4002, 4196, 4498, 4562, 4626, 4820, 5072]) = -2*u11(8)*X11(12);
	M([2080, 2469, 3681, 3940, 4004, 4133, 4197, 4261, 4693, 4757, 4821, 4885, 5137]) = -2*u11(7)*X11(12);
	M([3315, 3704, 4331, 4590, 4654, 4783, 4847, 4911, 4953, 5017, 5081, 5145, 5202]) = -u11(8)*X11(10)+u11(7)*X11(11);

	Mr = rref(M);  % replace me with a MEX

	A = zeros(16);
	amcols = [81 80 79 78 77 76 75 74 73 72 70 68 67 66 65 64];
	A(1, 4) = 1;
	A(2, 8) = 1;
	A(3, :) = -Mr(65, amcols);
	A(4, :) = -Mr(64, amcols);
	A(5, 12) = 1;
	A(6, 16) = 1;
	A(7, :) = -Mr(61, amcols);
	A(8, :) = -Mr(59, amcols);
	A(9, :) = -Mr(58, amcols);
	A(10, :) = -Mr(55, amcols);
	A(11, :) = -Mr(52, amcols);
	A(12, :) = -Mr(49, amcols);
	A(13, :) = -Mr(48, amcols);
	A(14, :) = -Mr(45, amcols);
	A(15, :) = -Mr(42, amcols);
	A(16, :) = -Mr(40, amcols);

	% generated for action variable 'd'
	% action matrix monomials: 1 b c d w b^2 c*b d*b w*b c^2 w*c w*d w^2 b^3 c*b^2 d*b^2
	[V D] = eig(A);
	sol =  V([5, 4, 3, 2],:)./(ones(4, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		w = [];
		d = [];
		c = [];
		b = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		w = sol(1,I);
		d = sol(2,I);
		c = sol(3,I);
		b = sol(4,I);
	end
end