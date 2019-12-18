from   math              import *
from   random            import *
import matplotlib.pyplot as plt
import copy

A  = (0.5) / 2
r  = 5
P  = 0.95
e  = 0.01
K  = 100
s1 = log(1 - P)
s2 = pi - 0
s3 = log(1 - e / s2)
N  = floor(s1 / s3)

def main(argv, argc):
	if argc == 1:
		Start(0, pi, 100000000)
	else:
		print("ERROR")

# Fg = ф с горизонтальной
# Fv - ф с волнистой

def get_W_A(r):
	# Функция определяет вес эелементов (Коэффициенты альфа)
	# Создаём массив на эр элементов (В нашем случае эр рвно 3)
	m    = floor((r - 1) / 2)
	a    = [0] * r
	a[m] = uniform(0, 1)
	# это центральный элемент, в него ставим рандомное число от 0 до 1
	Lsum, vseAll = a[m], a[m]
	# это левый столбик и центральный столбик для эр равно 3
	for i in range(m - 1, 0, -1):
		# запускам цикл для заполнения оставшихся элементов массива
		# 0 0,7 0
		# 0,5*ранд(0, 0.3)
		a[i] = 0.5 * uniform(0, 1 - Lsum)
		# 0,15 0,7 0
		a[r - 1 - i] = a[i]
		# 0.15 0.7 0.15
		Lsum   += a[i]
		vseAll += 2 * a[i]
	a[0], a[r - 1], vseAll = 0.5 * (1 - a[m]), 0.5 * (1 - a[m]), vseAll + 1 - a[m]
	for i in range(r):
		a[i] /= vseAll
	return a


def get_F_G(Fv, a, r, k):
	# реализуем формуЛУ расчета ФИЛтр СИГНАл по методу ср.арифм.
	memor = 0
	ss = (r-1)/2
	m   = floor(ss)
	Fg  = []
	for i in range(m, k - m + 1):
		# fk в методе равно fi у нас
		for j in range(i - m, i + m + 1):
			memor += a[j + m - i] * Fv[j]
		Fg.append(memor)
		memor = 0
	return Fg


def findWD(Fg, Fv, r, k):
	# находим критерйи зашум(W) и критерий отличия(D)
	m     = floor((r - 1) / 2)
	w, d  = [], []
	W = D = 0
	for i in range(0, k - m - 1):
		D += abs(Fg[i] - Fv[i])
	D *= 1 / k
	for i in range(1, k - m - 1):
		W += abs(Fg[i] - Fg[i - 1])
	return [W, D]


def Obnova(Fv, r, k):
	a          = get_W_A(r)
	Fg         = get_F_G(copy.deepcopy(Fv), copy.deepcopy(a), r, k)
	receivedWD = findWD(copy.deepcopy(Fg), copy.deepcopy(Fv), r, k)
	return [receivedWD[0], receivedWD[1], a, Fg]


def Start(Minx, Maxx, inf):
	global A, r, P, e, K, N
	x     = [0] * (K + 1)
	Fk    = [0] * (K + 1)
	Fv    = [0] * (K + 1)
	noise = [0] * (K + 1)
	L     = 10
	ss    = (r - 1) / 2
	M     = floor(ss)
	print("%-4s%8s%11s%8s%19s" % ("k", "x", "Fk", "noise", "Fv"))
	for i in range(K + 1):
		x[i]     = Minx + i * (Maxx - Minx) / K
		Fk[i]    = sin(x[i]) + 2 * A
		noise[i] = uniform(-A, A)
		Fv[i]    = Fk[i] + noise[i]
		print("%-4d%8.3f%9.4f%11.4f%16.4f" % (i, x[i], Fk[i], noise[i], Fv[i]))
	print()
	plt.title("Омега-Дельта")
	MinRast = inf
	for l in range(L + 1):
		ll = l / L
		print("Лямбда", ll)
		arrRast = Obnova(copy.deepcopy(Fv), r, K)
		for i in range(0, N):
			Rast = abs(arrRast[0]) + abs(arrRast[1])
			if i == 0:
				MinRast = Rast
				memor = arrRast
			else:
				if Rast < MinRast:
					MinRast = Rast
					memor = arrRast
		print("J( ", memor[0], "; ", memor[1], ") = ",
		      ll * memor[0] + (1 - ll) * memor[1])
		print("Rast =  ", memor[0])
		plt.plot(memor[0], memor[1], linestyle=":", marker='o', label=ll)
		if l == 0:
			memorMin = memor
			dmin   = MinRast
			llMin  = l / L
		else:
			if MinRast < dmin:
				memorMin = memor
				dmin   = MinRast
				llMin  = l / L
	plt.legend(loc='upper left')
	plt.savefig('C:/Users/danil/Desktop/mo/w-del1_r5.png', format='png')
	plt.show()
	plt.plot(range(0, K + 1), Fv, color='g', label="f_real")
	plt.plot(range(0, K + 1), Fk, color='b', label="f_ideal")
	if r == 5:
		plt.plot(range(K - M - 1), memorMin[3], color='r', label="result")
	if r == 3:
		plt.plot(range(K - M), memorMin[3], color='r', label="result")
	plt.savefig('C:/Users/danil/Desktop/mo/plot1_r5.png', format='png')
	plt.legend(loc='upper left')
	plt.show()
	plt.title("Noise")
	plt.plot(Fk, noise, color='y')
	plt.savefig('C:/Users/danil/Desktop/mo/noise1_r5.png', format='png')
	plt.show()
	print("\n\n\n")
	print("Результат")
	print("J = ",      llMin * memorMin[0] + (1 - llMin) * memorMin[1])
	print("w = ",      memorMin[0])
	print("o = ",      memorMin[1])
	print("alfa = ",   memorMin[2])
	print("lambda = ", llMin)


main(1, 1)
