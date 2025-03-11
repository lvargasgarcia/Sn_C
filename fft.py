import torch
import cnine
from fourierTransform import FourierTransform
import sys


def analysis(instance, output):
	f = open(output, "w")
	n = instance.getN()
	f1 = instance.getFunction()
	ft = FourierTransform(n, f1, mode="YSR", truncating_order=None)
	f.write(str(ft))
	f.close()

if __name__ == '__main__':
	from argparse import ArgumentParser,RawDescriptionHelpFormatter,_StoreTrueAction,ArgumentDefaultsHelpFormatter,Action
	parser = ArgumentParser(description = "Permutation surrogates")
	parser.add_argument('--problem', type=str, help='Problem: smwtp, samples')
	parser.add_argument('--instance', type=str, help='instance file')
	parser.add_argument("--output", type=str, default=None, help="output file")
	args = parser.parse_args()

	if args.problem == 'smwtp':
		from smwtp import SMWTP
		instance = SMWTP(args.instance)
	elif args.problem == 'samples':
		from functionsamples import FunctionFromSamples
		instance = FunctionFromSamples(args.instance)
	else:
		raise ValueError(f'Unsupported problem: {args.problem}')

	analysis(instance, args.output)
