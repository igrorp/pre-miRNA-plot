

desctxt = '''
 python3+ premiRNA-plot.py [INPUTS] ... [OPTIONS] ...

-----------------------------------------------------------------------------------------

Pre-miRNA-plot is a program for generating multiple custom images of miRNA precursors based
on RNAfold and RNAplot. It allows you to highlight the miRNA location within the precursor 
and obtain general and practical information about your data, so you can filter it or use 
it in publications.

-----------------------------------------------------------------------------------------

The program accepts tab-separated text files containing possibly 4 columns:

	1) Sequence ID (optional): Some sort of annotation or ID
		information about the sequence (e.g. 'ath-miR-171' or 'seq1');
	2) Precursor sequence: The pre-miRNA sequence.
	3) miRNA1 sequence: One of the miRNAs sequences.
	4) miRNA2 sequence (optional): The other miRNA sequence.
	
	File format examples:
	
	| pre-miRNA |   miRNA1   |  miRNA2  |
	
	>>>>>>>>>>>>>>>>>>>>>> OR <<<<<<<<<<<<<<<<<<<<<<<
	
	|     ID    |  pre-miRNA |  miRNA1  |   miRNA2  |

There is no problem if you don't have both miRNAs sequences; 
you can inform just one and the program will work just fine. 

Checkout our github repository for more details on how to use
the program: https://github.com/igrorp/pre-miRNA-plot

Parameters:

-i, --input       (str ...)
	 
	Inform the names of/paths to the files that you want to use.

-a, --annot       (str) --> (T or F)
	
	Informs if you have some sort of sequence ID, such as a miRNA
	family annotation (e.g.'ath-miRNA-171', 'seq1'), necessarily 
	on the first column, so that the generated image files can be 
	named according to that ID.

	Default = F (False)

-s, --style       (int) --> (between 1 and 5)

	The style of created images. Check the repository to see how they
	look like.

	Default = 3

-c, --color       (str str) OR (int int int int int int)
	
	You can choose which colors to paint the miRNA sequence within
	the precursor. Always provide the 5p and 3p colors, respectively.
	You can choose the predefined color names blue, red, green, purple,
	pink, yellow, cyan, white, black and orange; or you can inform the 
	RGB codes of the colors you want.

	Ex: '-c blue green' for blue 5p and green 3p
	Ex: '-c 255 255 0 153 0 204' for yellow 5p and purple 3p
	Default = green and red

-t, --threads     (int)
	
	Choose the number of allowed processors (CPUs) to be used.
	
	Default = 1

-f, --formats     (str) --> (svg or pdf)

	Choose the output format of the images. Choose between PDF and SVG.

	Default = svg

-o, --outdir      (str)

	The output name of directory created containing all the generated data.

	Default = premirnaplot



---------------------------------------------------------------

If you have any comments, complaints, doubts or suggestions, please contact
our main responsible for this project at igorpaim8@gmail.com or create an issue in
out github repository https://github.com/igrorp/pre-miRNA-plot/issues.

Have a nice work and let's keep making science evolve!

'''


defcolors = {
	'blue':'#0099ff',
	'red':'#ff0000',
	'green':'#33cc33',
	'purple':'#9933ff',
	'pink':'#ff66cc',
	'yellow':'#ffff00',
	'cyan':'#00ffff',
	'white':'#ffffff',
	'black':'#000000',
	'orange':'#ff9933'
}

# RNAfold precalculated values of MFE expected for RNA sequences with L nucleotides and equimolar ratios of A, C, G and U 

refmfe = {
	40:	-6.9620,
	41:	-7.2395,
	42:	-7.5170,
	43:	-7.7945,
	44:	-8.0720,
	45:	-8.3305,
	46:	-8.5890,
	47:	-8.8475,
	48:	-9.1060,
	49:	-9.3445,
	50:	-9.5830,
	51:	-9.8215,
	52:	-10.0600,
	53:	-10.3833,
	54:	-10.7065,
	55:	-11.0298,
	56:	-11.3530,
	57:	-11.6108,
	58:	-11.8685,
	59:	-12.1263,
	60:	-12.3840,
	61:	-12.6833,
	62:	-12.9825,
	63:	-13.2818,
	64:	-13.5810,
	65:	-13.8590,
	66:	-14.1370,
	67:	-14.4150,
	68:	-14.6930,
	69:	-14.9810,
	70:	-15.2690,
	71:	-15.5570,
	72:	-15.8450,
	73:	-16.1090,
	74:	-16.3730,
	75:	-16.6370,
	76:	-16.9010,
	77:	-17.1763,
	78:	-17.4515,
	79:	-17.7268,
	80:	-18.0020,
	81:	-18.2925,
	82:	-18.5830,
	83:	-18.8735,
	84:	-19.1640,
	85:	-19.4395,
	86:	-19.7150,
	87:	-19.9905,
	88:	-20.2660,
	89:	-20.5860,
	90:	-20.9060,
	91:	-21.2260,
	92:	-21.5460,
	93:	-21.8253,
	94:	-22.1045,
	95:	-22.3838,
	96:	-22.6630,
	97:	-22.9770,
	98:	-23.2910,
	99:	-23.6050,
	100:	-23.9190,
	101:	-24.1945,
	102:	-24.4700,
	103:	-24.7455,
	104:	-25.0210,
	105:	-25.3165,
	106:	-25.6120,
	107:	-25.9075,
	108:	-26.2030,
	109:	-26.4780,
	110:	-26.7530,
	111:	-27.0280,
	112:	-27.3030,
	113:	-27.5890,
	114:	-27.8750,
	115:	-28.1610,
	116:	-28.4470,
	117:	-28.7663,
	118:	-29.0855,
	119:	-29.4048,
	120:	-29.7240,
	121:	-30.0300,
	122:	-30.3360,
	123:	-30.6420,
	124:	-30.9480,
	125:	-31.2528,
	126:	-31.5575,
	127:	-31.8623,
	128:	-32.1670,
	129:	-32.4330,
	130:	-32.6990,
	131:	-32.9650,
	132:	-33.2310,
	133:	-33.5518,
	134:	-33.8725,
	135:	-34.1933,
	136:	-34.5140,
	137:	-34.8333,
	138:	-35.1525,
	139:	-35.4718,
	140:	-35.7910,
	141:	-36.1063,
	142:	-36.4215,
	143:	-36.7368,
	144:	-37.0520,
	145:	-37.2758,
	146:	-37.4995,
	147:	-37.7233,
	148:	-37.9470,
	149:	-38.3218,
	150:	-38.6965,
	151:	-39.0713,
	152:	-39.4460,
	153:	-39.7530,
	154:	-40.0600,
	155:	-40.3670,
	156:	-40.6740,
	157:	-41.0090,
	158:	-41.3440,
	159:	-41.6790,
	160:	-42.0140,
	161:	-42.2888,
	162:	-42.5635,
	163:	-42.8383,
	164:	-43.1130,
	165:	-43.4003,
	166:	-43.6875,
	167:	-43.9748,
	168:	-44.2620,
	169:	-44.5498,
	170:	-44.8375,
	171:	-45.1253,
	172:	-45.4130,
	173:	-45.7080,
	174:	-46.0030,
	175:	-46.2980,
	176:	-46.5930,
	177:	-46.8978,
	178:	-47.2025,
	179:	-47.5073,
	180:	-47.8120,
	181:	-48.1605,
	182:	-48.5090,
	183:	-48.8575,
	184:	-49.2060,
	185:	-49.5183,
	186:	-49.8305,
	187:	-50.1428,
	188:	-50.4550,
	189:	-50.7805,
	190:	-51.1060,
	191:	-51.4315,
	192:	-51.7570,
	193:	-52.0605,
	194:	-52.3640,
	195:	-52.6675,
	196:	-52.9710,
	197:	-53.2478,
	198:	-53.5245,
	199:	-53.8013,
	200:	-54.0780,
	201:	-54.3628,
	202:	-54.6475,
	203:	-54.9323,
	204:	-55.2170,
	205:	-55.5388,
	206:	-55.8605,
	207:	-56.1823,
	208:	-56.5040,
	209:	-56.8073,
	210:	-57.1105,
	211:	-57.4138,
	212:	-57.7170,
	213:	-58.0323,
	214:	-58.3475,
	215:	-58.6628,
	216:	-58.9780,
	217:	-59.3283,
	218:	-59.6785,
	219:	-60.0288,
	220:	-60.3790,
	221:	-60.6730,
	222:	-60.9670,
	223:	-61.2610,
	224:	-61.5550,
	225:	-61.8765,
	226:	-62.1980,
	227:	-62.5195,
	228:	-62.8410,
	229:	-63.0798,
	230:	-63.3185,
	231:	-63.5573,
	232:	-63.7960,
	233:	-64.1090,
	234:	-64.4220,
	235:	-64.7350,
	236:	-65.0480,
	237:	-65.3523,
	238:	-65.6565,
	239:	-65.9608,
	240:	-66.2650,
	241:	-66.6333,
	242:	-67.0015,
	243:	-67.3698,
	244:	-67.7380,
	245:	-68.0510,
	246:	-68.3640,
	247:	-68.6770,
	248:	-68.9900,
	249:	-69.3105,
	250:	-69.6310,
	251:	-69.9515,
	252:	-70.2720,
	253:	-70.5268,
	254:	-70.7815,
	255:	-71.0363,
	256:	-71.2910,
	257:	-71.6643,
	258:	-72.0375,
	259:	-72.4108,
	260:	-72.7840,
	261:	-73.0830,
	262:	-73.3820,
	263:	-73.6810,
	264:	-73.9800,
	265:	-74.2858,
	266:	-74.5915,
	267:	-74.8973,
	268:	-75.2030,
	269:	-75.5108,
	270:	-75.8185,
	271:	-76.1263,
	272:	-76.4340,
	273:	-76.7523,
	274:	-77.0705,
	275:	-77.3888,
	276:	-77.7070,
	277:	-78.0683,
	278:	-78.4295,
	279:	-78.7908,
	280:	-79.1520,
	281:	-79.4495,
	282:	-79.7470,
	283:	-80.0445,
	284:	-80.3420,
	285:	-80.5848,
	286:	-80.8275,
	287:	-81.0703,
	288:	-81.3130,
	289:	-81.5858,
	290:	-81.8585,
	291:	-82.1313,
	292:	-82.4040,
	293:	-82.7833,
	294:	-83.1625,
	295:	-83.5418,
	296:	-83.9210,
	297:	-84.2398,
	298:	-84.5585,
	299:	-84.8773,
	300:	-85.1960,
	301:	-85.4578,
	302:	-85.7195,
	303:	-85.9813,
	304:	-86.2430,
	305:	-86.6260,
	306:	-87.0090,
	307:	-87.3920,
	308:	-87.7750,
	309:	-88.0850,
	310:	-88.3950,
	311:	-88.7050,
	312:	-89.0150,
	313:	-89.3113,
	314:	-89.6075,
	315:	-89.9038,
	316:	-90.2000,
	317:	-90.4365,
	318:	-90.6730,
	319:	-90.9095,
	320:	-91.1460,
	321:	-91.5115,
	322:	-91.8770,
	323:	-92.2425,
	324:	-92.6080,
	325:	-92.9770,
	326:	-93.3460,
	327:	-93.7150,
	328:	-94.0840,
	329:	-94.3583,
	330:	-94.6325,
	331:	-94.9068,
	332:	-95.1810,
	333:	-95.5083,
	334:	-95.8355,
	335:	-96.1628,
	336:	-96.4900,
	337:	-96.8243,
	338:	-97.1585,
	339:	-97.4928,
	340:	-97.8270,
	341:	-98.1080,
	342:	-98.3890,
	343:	-98.6700,
	344:	-98.9510,
	345:	-99.2468,
	346:	-99.5425,
	347:	-99.8383,
	348:	-100.1340,
	349:	-100.4603,
	350:	-100.7865,
	351:	-101.1128,
	352:	-101.4390,
	353:	-101.7790,
	354:	-102.1190,
	355:	-102.4590,
	356:	-102.7990,
	357:	-103.0393,
	358:	-103.2795,
	359:	-103.5198,
	360:	-103.7600,
	361:	-104.1548,
	362:	-104.5495,
	363:	-104.9443,
	364:	-105.3390,
	365:	-105.6178,
	366:	-105.8965,
	367:	-106.1753,
	368:	-106.4540,
	369:	-106.7870,
	370:	-107.1200,
	371:	-107.4530,
	372:	-107.7860,
	373:	-108.0905,
	374:	-108.3950,
	375:	-108.6995,
	376:	-109.0040,
	377:	-109.2688,
	378:	-109.5335,
	379:	-109.7983,
	380:	-110.0630,
	381:	-110.4093,
	382:	-110.7555,
	383:	-111.1018,
	384:	-111.4480,
	385:	-111.7515,
	386:	-112.0550,
	387:	-112.3585,
	388:	-112.6620,
	389:	-112.9708,
	390:	-113.2795,
	391:	-113.5883,
	392:	-113.8970,
	393:	-114.3525,
	394:	-114.8080,
	395:	-115.2635,
	396:	-115.7190,
	397:	-115.8665,
	398:	-116.0140,
	399:	-116.1615,
	400:	-116.3090,
	401:	-116.7353,
	402:	-117.1615,
	403:	-117.5878,
	404:	-118.0140,
	405:	-118.3453,
	406:	-118.6765,
	407:	-119.0078,
	408:	-119.3390,
	409:	-119.5980,
	410:	-119.8570,
	411:	-120.1160,
	412:	-120.3750,
	413:	-120.7270,
	414:	-121.0790,
	415:	-121.4310,
	416:	-121.7830,
	417:	-122.0363,
	418:	-122.2895,
	419:	-122.5428,
	420:	-122.7960,
	421:	-123.1490,
	422:	-123.5020,
	423:	-123.8550,
	424:	-124.2080,
	425:	-124.5475,
	426:	-124.8870,
	427:	-125.2265,
	428:	-125.5660,
	429:	-125.8953,
	430:	-126.2245,
	431:	-126.5538,
	432:	-126.8830,
	433:	-127.1908,
	434:	-127.4985,
	435:	-127.8063,
	436:	-128.1140,
	437:	-128.3733,
	438:	-128.6325,
	439:	-128.8918,
	440:	-129.1510,
	441:	-129.5333,
	442:	-129.9155,
	443:	-130.2978,
	444:	-130.6800,
	445:	-130.9718,
	446:	-131.2635,
	447:	-131.5553,
	448:	-131.8470,
	449:	-132.2133,
	450:	-132.5795,
	451:	-132.9458,
	452:	-133.3120,
	453:	-133.6295,
	454:	-133.9470,
	455:	-134.2645,
	456:	-134.5820,
	457:	-134.8568,
	458:	-135.1315,
	459:	-135.4063,
	460:	-135.6810,
	461:	-135.9910,
	462:	-136.3010,
	463:	-136.6110,
	464:	-136.9210,
	465:	-137.2958,
	466:	-137.6705,
	467:	-138.0453,
	468:	-138.4200,
	469:	-138.6513,
	470:	-138.8825,
	471:	-139.1138,
	472:	-139.3450,
	473:	-139.7158,
	474:	-140.0865,
	475:	-140.4573,
	476:	-140.8280,
	477:	-141.0943,
	478:	-141.3605,
	479:	-141.6268,
	480:	-141.8930,
	481:	-142.1658,
	482:	-142.4385,
	483:	-142.7113,
	484:	-142.9840,
	485:	-143.3620,
	486:	-143.7400,
	487:	-144.1180,
	488:	-144.4960,
	489:	-144.8195,
	490:	-145.1430,
	491:	-145.4665,
	492:	-145.7900,
	493:	-146.0723,
	494:	-146.3545,
	495:	-146.6368,
	496:	-146.9190,
	497:	-147.2530,
	498:	-147.5870,
	499:	-147.9210,
	500:	-148.2550,
	501:	-148.6025,
	502:	-148.9500,
	503:	-149.2975,
	504:	-149.6450,
	505:	-149.9278,
	506:	-150.2105,
	507:	-150.4933,
	508:	-150.7760,
	509:	-151.0550,
	510:	-151.3340,
	511:	-151.6130,
	512:	-151.8920,
	513:	-152.2803,
	514:	-152.6685,
	515:	-153.0568,
	516:	-153.4450,
	517:	-153.7778,
	518:	-154.1105,
	519:	-154.4433,
	520:	-154.7760,
	521:	-155.0910,
	522:	-155.4060,
	523:	-155.7210,
	524:	-156.0360,
	525:	-156.3630,
	526:	-156.6900,
	527:	-157.0170,
	528:	-157.3440,
	529:	-157.6468,
	530:	-157.9495,
	531:	-158.2523,
	532:	-158.5550,
	533:	-158.8875,
	534:	-159.2200,
	535:	-159.5525,
	536:	-159.8850,
	537:	-160.2405,
	538:	-160.5960,
	539:	-160.9515,
	540:	-161.3070,
	541:	-161.5970,
	542:	-161.8870,
	543:	-162.1770,
	544:	-162.4670,
	545:	-162.7640,
	546:	-163.0610,
	547:	-163.3580,
	548:	-163.6550,
	549:	-163.9650,
	550:	-164.2750,
	551:	-164.5850,
	552:	-164.8950,
	553:	-165.1808,
	554:	-165.4665,
	555:	-165.7523,
	556:	-166.0380,
	557:	-166.3568,
	558:	-166.6755,
	559:	-166.9943,
	560:	-167.3130,
	561:	-167.6345,
	562:	-167.9560,
	563:	-168.2775,
	564:	-168.5990,
	565:	-168.9230,
	566:	-169.2470,
	567:	-169.5710,
	568:	-169.8950,
	569:	-170.2070,
	570:	-170.5190,
	571:	-170.8310,
	572:	-171.1430,
	573:	-171.4415,
	574:	-171.7400,
	575:	-172.0385,
	576:	-172.3370,
	577:	-172.7660,
	578:	-173.1950,
	579:	-173.6240,
	580:	-174.0530,
	581:	-174.2958,
	582:	-174.5385,
	583:	-174.7813,
	584:	-175.0240,
	585:	-175.3080,
	586:	-175.5920,
	587:	-175.8760,
	588:	-176.1600,
	589:	-176.5448,
	590:	-176.9295,
	591:	-177.3143,
	592:	-177.6990,
	593:	-178.0120,
	594:	-178.3250,
	595:	-178.6380,
	596:	-178.9510,
	597:	-179.2725,
	598:	-179.5940,
	599:	-179.9155,
	600:	-180.2370,
}

# Predefined optimal constant amount that shifts the MFE x length regression line to the origin of the graph

SHIFT_CONST = 8


class Precursor():

	def __init__(self, name, mirna1, mirna2, precursor):

		# The pre,fix of the filename for the pre-miRNA
		
		self.name = name

		# The precursor nucleotide sequence
		
		self.premirna = precursor

		# A tuple containing the miRNAs
		
		self.mirnas = (mirna1 if mirna1 else '', mirna2 if mirna2 else '')
		
		# The tuple containing the positions of the miRNA within the pre-miRNA or None
		
		self.pos1 = self.__pos(mirna1)

		# The tuple containing the positions of the other miRNA within the pre-miRNA or None
		
		self.pos2 = self.__pos(mirna2)

		# The length of the pre-miRNA
		
		self.prelen = len(precursor)

		# The Minimum Free Energy for the precursor as predicted by RNAfold
		
		self.premfe = 0

		# The precursor's predict secondary structure in dot bracket notation
		
		self.predsec = ''

		# The GC content of the precursor

		self.gccontent = self.__gc()

		# The MFEdensity of the precursor

		self.mfeden = self.__mfeden()

		# The number of mismatches in the region of the miRNA duplex

		self.mm = 'N.A.'

		# The number of bulges from the precursor


	def __pos(self, mirna):

		if mirna:

			if mirna in self.premirna:

				if self.premirna.count(mirna) > 1:
					
					print('WARNING! miRNA {} was found more than once in the precursor sequence {}..., but its last occurrence will be used!'.format(mirna, self.premirna[25:]))

				return (self.premirna.find(mirna), self.premirna.find(mirna) + len(mirna))
			
			else:
				
				raise Exception('ERROR! Could not find sequence {} inside {}, please correct this'.format(mirna, self.premirna))

		else:

			return None

	
	def __mfeden(self):

		if self.prelen >= 40 and self.prelen <= 600:

			return round(100 * (self.premfe - refmfe[self.prelen]) / (self.prelen - SHIFT_CONST), 2)

		else:

			return 'N.A.'


	def __gc(self):

		return round((self.premirna.count('G') + self.premirna.count('C')) / self.prelen , 2)


	def mismatches(self):

		if self.pos1 and self.pos2:

			posa, posb = self.pos1

			posc, posd = self.pos2

			self.mm = self.predsec[posa:posb].count('.') + self.predsec[posc:posd].count('.')

		else:

			return 'N.A.'

	def setpremfe(self, mfe):

		self.premfe = mfe

		self.mfeden = self.__mfeden()



		


		