import math
import itertools
import numpy as np
import fourier_transform as fourier
from smwtp import *
from fractions import Fraction

global TOLERANCE

TOLERANCE = 0.001

def maeMaxMin(f1, f2, n):
	sum, f1Min, f1Max, f2Min, f2Max = (0,None,None,None,None)
	for permutation in itertools.permutations(range(1,n+1)):
		valF1 = f1(permutation)
		valF2 = f2(permutation)
		if f1Min == None or valF1 < f1Min:
			f1Min = valF1
		if f1Max == None or valF1 > f1Max:
			f1Max = valF1
		if f2Min == None or valF2 < f2Min:
			f2Min = valF2
		if f2Max == None or valF2 > f2Max:
			f2Max = valF2
		sum = sum + abs(valF1-valF2)
	return (sum/math.factorial(n), f1Min, f1Max, f2Min, f2Max)

def maeOfGlobalOptima(f1, f2, n):
	sumF1 = 0
	globalOptimaValueF1 = None
	globalOptimaF1 = 0
	sumF2 = 0
	globalOptimaValueF2 = None
	globalOptimaF2 = 0

	globalOptimaF2List = list()

	for permutation in itertools.permutations(range(1,n+1)):
		valF1 = f1(permutation)
		valF2 = f2(permutation)
		if globalOptimaValueF1 == None or valF1 < globalOptimaValueF1:
			globalOptimaValueF1 = valF1
			sumF1 = abs(valF1-valF2)
			globalOptimaF1 = 1
		elif valF1 == globalOptimaValueF1:
			sumF1 = sumF1 + abs(valF1-valF2)
			globalOptimaF1 = globalOptimaF1 + 1

		if globalOptimaValueF2 == None or valF2 < globalOptimaValueF2:
			globalOptimaValueF2 = valF2
			sumF2 = abs(valF1-valF2)
			globalOptimaF2 = 1
			globalOptimaF2List = [permutation]
		elif valF2 == globalOptimaValueF2:
			sumF2 = sumF2 + abs(valF1-valF2)
			globalOptimaF2 = globalOptimaF2 + 1
			globalOptimaF2List.append(permutation)

	preservedGlobalOptima = len([permutation for permutation in globalOptimaF2List if f1(permutation)==globalOptimaValueF1])


	return sumF1/globalOptimaF1, globalOptimaF1, sumF2/globalOptimaF2, globalOptimaF2, preservedGlobalOptima

def allPermutations(n):
	return [permutation for permutation in itertools.permutations(range(1, n + 1))]

def sortedPermutations(n, f):
	list = allPermutations(n)
	list.sort(key=(lambda p: f(p)))
	return list

def permutationRanking(n, f):
	permutations = sortedPermutations(n,f)
	result = {}

	previousValue = None
	rank = 0
	for p in permutations:
		if previousValue == None or f(p) > previousValue + TOLERANCE:
			rank = rank+1
			previousValue = f(p)
		result[tuple(p)] = rank
	
	def resp(p):
		return result[p]

	return resp


def showNormalOrder(n, f):
	print("Normal order")
	for p in allPermutations(n):
		print(f'{p}\t{f[p]}')


def showSortedOrder(n, f):
	print("Sorted permutations")
	for p in sortedPermutations(n, f):
		print(f'{p}\t{f[p]}')


def showPermutationRanking(n, f):
	print("Permutation ranking")
	ranking = permutationRanking(n, f)
	for p in allPermutations(n):
		print(f'{p}\t{ranking[p]}')

def analysis(instance, output, mode="YSR", type="double"):
	
	if type == "double":
		f = open(output, "w")
		f.write("----------------- Floating point (64 bits) -----------------\n")
	else:
		f = open(output, "a")
		f.write("----------------- Rationals -------------------\n")

	n = instance.getN()
	dict_instance = instance.getFunction()

	if type == "rational":
		for p in dict_instance:
			dict_instance[p] = fourier.Rational(int(dict_instance[p]),1)

	if type == "double":
		ft = fourier.FourierTransform_double(n, dict_instance, mode, 6)
	else:
		ft = fourier.FourierTransform_Rational(n, dict_instance, mode, 6)
	
	def f1(p):
		fr = dict_instance[to_int(p)]
		return fr if isinstance(fr,float) else Fraction(fr.numerator, fr.denominator)

	#showNormalOrder(n, f1)
	#showSortedOrder(n, f1)
	#showPermutationRanking(n, f1)
	#print(f'Global min: {instance.globalMin}')
	#print(f'Global max: {instance.globalMax}')
	#print("f1")
	#showSortedOrder(n, f1)
	rankingF1 = permutationRanking(n, f1)
	#print("Ranking f1")
	#showNormalOrder(n, rankingF1)
	f.write('Max order\tMAE\tNormalized MAE\tF1 Min\tF1 Max\tF2 Min\tF2 Max\tMAE-GO Orig\tNormalized MAE-GO Orig\tGO Orig\tMAE-GO Trunc\tNormalized MAE-GO Trunc\tGO Trunc\tRanking MAE\tPreserved GO\n')
	for firstLine in range(0,n):
		
		def f2(p):
			fr = ft.inverseFT(p, firstLine)
			return fr if isinstance(fr, float) else Fraction(fr.numerator,fr.denominator)
		
		#showPermutationRanking(n, f2)
		val, f1Min, f1Max, f2Min, f2Max = maeMaxMin(f1, f2, n)
		f2Min = (float)(f2Min)
		f2Max = (float)(f2Max)
		fRange = (instance.globalMax-instance.globalMin)
		maeGOF1, globalOptimaF1, maeGOF2, globalOptimaF2, preservedGlobalOptima = maeOfGlobalOptima(f1, f2, n)

		rankingF2 = permutationRanking(n, f2)
		#print("f2")
		#showSortedOrder(n, f2)
		#print("Ranking f2")
		#showSortedOrder(n,rankingF2)
		maeRanking, _, _, _, _ = maeMaxMin(rankingF1, rankingF2, n)
		f.write(f'{firstLine}\t{val}\t{val/fRange}\t{f1Min}\t{f1Max}\t{f2Min}\t{f2Max}\t{maeGOF1}\t{maeGOF1 / fRange}\t{globalOptimaF1}\t{maeGOF2}\t{maeGOF2 / fRange}\t{globalOptimaF2}\t{maeRanking}\t{preservedGlobalOptima}\n')
	
	f.close()

if __name__ == '__main__':
	from argparse import ArgumentParser,RawDescriptionHelpFormatter,_StoreTrueAction,ArgumentDefaultsHelpFormatter,Action
	parser = ArgumentParser(description = "Permutation surrogates")
	parser.add_argument('--problem', type=str, help='Problem: smwtp, samples')
	parser.add_argument('--instance', type=str, help='instance file')
	parser.add_argument("--output", type=str, default=None, help="output file")
	parser.add_argument("--mode", type=str, default="YSR", help="mode")
	args = parser.parse_args()

	if args.problem == 'smwtp':
		from smwtp import SMWTP
		instance = SMWTP(args.instance)
	elif args.problem == 'samples':
		from functionsamples import FunctionFromSamples
		instance = FunctionFromSamples(args.instance)
	else:
		raise ValueError(f'Unsupported problem: {args.problem}')
	
	analysis(instance, args.output, args.mode, "double")
	analysis(instance, args.output, args.mode, "rational")
