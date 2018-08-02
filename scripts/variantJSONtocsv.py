import argparse
import json
import csv

def parseJson(data):
	parsedDataSet = []
	sample = data["Sample"]
	genome = data["genome"]
	for i in range(0,len(genome)):
		#cycle through each segment
		chr = genome[i]["chr"]
		for j in range(0,len(genome[i]['seq'])):
			consensus = genome[i]['seq'][j]['consensus']
			coverage = genome[i]['seq'][j]['coverage']
			concat_pos = genome[i]['seq'][j]['concat_pos']
			pos = genome[i]['seq'][j]['pos']
			for key in genome[i]['seq'][j]['alleles'].keys():
				freq = genome[i]['seq'][j]['alleles'][key]['freq']
				nt = genome[i]['seq'][j]['alleles'][key]['nucleotide']
				count = genome[i]['seq'][j]['alleles'][key]['count']
				datapoint = {
					"Sample": sample,
					"chr": chr,
					"nucleotide": nt,
					"consensus": consensus,
					"pos": pos,
					"concat_pos": concat_pos,
					"freq": freq,
					"count": count,
					"coverage": coverage,
					"mutationalClass": genome[i]['seq'][j]['alleles'][key]['mutationalClass']
				}
				parsedDataSet.append(datapoint)
		return(parsedDataSet)


def main():
	parser = argparse.ArgumentParser(description='This scipts takes a json file as outputed by the pipeline and converts it to a csv file with one allele per line',
	usage ="python jsonToCsv.py sample.json sample.csv")
	parser.add_argument('json_file', meta='json_file', nargs='+',
                        help='a json file from the pipeline')

	parser.add_argument('output_csv', meta='output_csv',nargs='+',
                            help = 'the output csv')

	args = parser.parse_args()

	with open(args.json_file[0]) as json_file:
		data = json.load(json_file)

	parsedData = parseJson(data)
	headers=parsedData[0].keys()
	with open(args.output_csv[0], 'wb') as output_file:
		dict_writer = csv.DictWriter(output_file, headers)
		dict_writer.writeheader()
		dict_writer.writerows(parsedData)
	



if __name__ == '__main__':
	main()


