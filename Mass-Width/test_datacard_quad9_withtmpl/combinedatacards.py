

import argparse
quad = 9

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-o','--output',dest='outputcard',default='hzz4l_combined_13TeV_onshell.txt')
	parser.add_argument('-c','--channels',action='append',default=['all'])
	parser.add_argument('--cat',action='append',default=['all'])
	args = parser.parse_args()
	
	#hzz4l_4e_cat1_13TeV_onshell.txt
	if 'all' in args.channels:		channels = ['4e','4mu','2e2mu']
	else:					channels = args.channels
	if 'all' in args.cat:			cats = [str(cat) for cat in range(quad)]
	else:					cats = args.cat
	
	combcmd = 'combineCards.py '
	for ch in channels:
		for cat in cats:
			combcmd += 'hzz4l_{ch}_cat{cat}_13TeV_onshell.txt '.format(ch=ch,cat=cat)
	combcmd += ' > {output}'.format(output=args.outputcard)
	os.system(combcmd)
