import math
import itertools
import numpy as np
import fourier_transform as fourier
from smwtp import *
from fractions import Fraction
import json
import os 

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

def maeOfGlobalOptima(f1, f2, n, epsilon):
	sumF1 = 0
	globalOptimaValueF1 = None
	globalOptimaF1 = 0
	sumF2 = 0
	globalOptimaValueF2 = None
	globalOptimaF2 = 0

	globalOptimaF1List = list()
	globalOptimaF2List = list()

	for permutation in itertools.permutations(range(1,n+1)):
		valF1 = f1(permutation)
		valF2 = f2(permutation)
		if globalOptimaValueF1 == None or valF1 < globalOptimaValueF1:
			globalOptimaValueF1 = valF1
			globalOptimaF1List = [permutation]
			sumF1 = abs(valF1-valF2)
			globalOptimaF1 = 1
		elif valF1 == globalOptimaValueF1:
			sumF1 = sumF1 + abs(valF1-valF2)
			globalOptimaF1 = globalOptimaF1 + 1
			globalOptimaF1List.append(permutation)

		if globalOptimaValueF2 == None or globalOptimaValueF2 - valF2 > epsilon:
			globalOptimaValueF2 = valF2
			sumF2 = abs(valF1-valF2)
			globalOptimaF2 = 1
			globalOptimaF2List = [permutation]
		elif abs(valF2 - globalOptimaValueF2) <= epsilon:
			sumF2 = sumF2 + abs(valF1-valF2)
			globalOptimaF2 = globalOptimaF2 + 1
			globalOptimaF2List.append(permutation)

	# print(globalOptimaF1List)
	# print(globalOptimaF2List)
	# print("-----------------------------")
	preservedGlobalOptima = len([permutation for permutation in globalOptimaF2List if f1(permutation)==globalOptimaValueF1])


	return sumF1/globalOptimaF1, globalOptimaF1, sumF2/globalOptimaF2, globalOptimaF2, preservedGlobalOptima

def distance_of_GOs(f1,f2,n):
	globalOptimaValueF1 = None
	globalOptimaValueF2 = None
	optimal_pi_F2 = None
	for pi in itertools.permutations(range(1,n+1)):
		valF1 = f1(pi)
		valF2 = f2(pi)
		if globalOptimaValueF1 == None or valF1 < globalOptimaValueF1:
			globalOptimaValueF1 = valF1
		if globalOptimaValueF2 == None or valF2 < globalOptimaValueF2:
			globalOptimaValueF2 = valF2
			optimal_pi_F2 = pi
	return abs(globalOptimaValueF1-f1(optimal_pi_F2))

def allPermutations(n):
	return [permutation for permutation in itertools.permutations(range(1, n + 1))]

def sortedPermutations(n, f):
	list = allPermutations(n)
	list.sort(key=(lambda p: f(p)))
	return list

def permutationRanking(n, f, epsilon):
	permutations = sortedPermutations(n,f)
	result = {}

	previousValue = None
	rank = 0
	for p in permutations:
		if previousValue == None or f(p) > previousValue + epsilon:
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

def analysis(epsilon, name, instance, mode="YSR", type="double"):
	
	info = {}

	info["name"] = name + ":" + type
	info["mode"] = mode
	info["tolerance"] = epsilon

	info["orders"] = {}

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
		return fr if (isinstance(fr,float) or isinstance(fr, int)) else Fraction(fr.numerator, fr.denominator)

	rankingF1 = permutationRanking(n, f1, epsilon)

	for firstLine in range(0,n):
		
		info["orders"][firstLine] = {}

		def f2(p):
			fr = ft.inverseFT(p, firstLine)
			return fr if isinstance(fr, float) else Fraction(fr.numerator,fr.denominator)
		
		#showPermutationRanking(n, f2)
		val, f1Min, f1Max, f2Min, f2Max = maeMaxMin(f1, f2, n)
		f2Min = (float)(f2Min)
		f2Max = (float)(f2Max)
		fRange = (instance.globalMax-instance.globalMin)
		maeGOF1, globalOptimaF1, maeGOF2, globalOptimaF2, preservedGlobalOptima = maeOfGlobalOptima(f1, f2, n, epsilon)

		rankingF2 = permutationRanking(n, f2, epsilon)

		maeRanking, _, _, _, _ = maeMaxMin(rankingF1, rankingF2, n)
		distance_gos = float(distance_of_GOs(f1,f2,n))
		val = float(val)
		maeGOF1 = float(maeGOF1)
		maeGOF2 = float(maeGOF2)
		globalOptimaF1 = float(globalOptimaF1)
		globalOptimaF2 = float(globalOptimaF2)
		maeRanking = float(maeRanking)
		info["orders"][firstLine]["MAE"] = val
		info["orders"][firstLine]["Normalized MAE"] = val/fRange
		info["orders"][firstLine]["F1 Min"] = f1Min
		info["orders"][firstLine]["F1 Max"] = f1Max
		info["orders"][firstLine]["F2 Min"] = f2Min
		info["orders"][firstLine]["F2 Max"] = f2Max
		info["orders"][firstLine]["MAE-GO Orig"] = maeGOF1
		info["orders"][firstLine]["Normalized MAE-GO Orig"] = maeGOF1 / fRange
		info["orders"][firstLine]["GO Orig"] = globalOptimaF1
		info["orders"][firstLine]["MAE-GO Trunc"] = maeGOF2
		info["orders"][firstLine]["Normalized MAE-GO Trunc"] = maeGOF2 / fRange
		info["orders"][firstLine]["GO Trunc"] = globalOptimaF2
		info["orders"][firstLine]["Ranking MAE"] = maeRanking
		info["orders"][firstLine]["Preserved GO"] = preservedGlobalOptima
		info["orders"][firstLine]["NPGO"] = preservedGlobalOptima/globalOptimaF1
		info["orders"][firstLine]["Distance GOs"] = distance_gos
		info["orders"][firstLine]["Normalized distance GOs"] = distance_gos/fRange

		# f.write(f'{firstLine}\t{val}\t{val/fRange}\t{f1Min}\t{f1Max}\t{f2Min}\t{f2Max}\t{maeGOF1}\t{maeGOF1 / fRange}\t{globalOptimaF1}\t{maeGOF2}\t{maeGOF2 / fRange}\t{globalOptimaF2}\t{maeRanking}\t{preservedGlobalOptima}\n')
	
	#f.close()
	output_file = "./results/" + mode + "/" + ("smwtp" if isinstance(instance, SMWTP) else "arp") + "/" + str(n) + "/"
	name = os.path.splitext(name.split("/")[-1])[0]
	json.dump(info, open(output_file + name + ("_fp" if type == "double" else "_r") + "_" + str(epsilon) + ".json", "w"), indent=4, default=str)

if __name__ == '__main__':
	from argparse import ArgumentParser,RawDescriptionHelpFormatter,_StoreTrueAction,ArgumentDefaultsHelpFormatter,Action
	parser = ArgumentParser(description = "Permutation surrogates")
	parser.add_argument('--problem', type=str, help='Problem: smwtp, samples')
	parser.add_argument('--instance', type=str, help='instance file')
	parser.add_argument("--mode", type=str, default="YSR", help="mode")
	args = parser.parse_args()

	if args.problem == 'smwtp':
		from smwtp import SMWTP
		instance = SMWTP(args.instance)
	elif args.problem == 'arp':
		from functionsamples import FunctionFromSamples
		instance = FunctionFromSamples(args.instance)
	else:
		raise ValueError(f'Unsupported problem: {args.problem}')

	# analysis(0.5, "instances/arp/arp_7_17.csv", FunctionFromSamples("instances/arp/arp_7_17.csv"), "r", "YSR", "double")
	
	analysis(0,args.instance, instance, args.mode, "double")
	analysis(0.5,args.instance, instance, args.mode, "double")
	if args.mode != "YOR":
		analysis(0,args.instance, instance, args.mode, "rational")
