#!/usr/bin/env python3
"""
Get eROSITA skytile. Note, this can be significantly slower than the local method, depending on your network.

"""
import sys
import argparse
from erotools.radec2tile import tile_api

class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2024 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('ra',type=str,help='source ra')
parser.add_argument('dec',type=str,help='source dec')
parser.add_argument('--radius',type=str,default=0,help='match radius')

args = parser.parse_args()


def main():
    skytile = tile_api(args.ra,args.dec,args.radius)
    print(skytile)

# Example usage
if __name__ == "__main__":
    main()
    

