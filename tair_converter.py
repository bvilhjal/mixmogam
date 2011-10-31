import bisect
import env

class tair8_to_tair9_map():
	"""
	Converts TAIR 8 positions on the given chromosome to the corresponding TAIR 9 positions. 
	"""
	def __init__(self):
		#Read conversion manual.
		file_name = env.env['data_dir'] + 'tair/TAIR9_assembly_updates_relative_to_TAIR8_TAIR9_assemblies.csv'
		self.offsets = {}
		self.offset_positions = {}
		for chrom in [1, 2, 3, 4, 5]:
			self.offsets[chrom] = []
			self.offset_positions[chrom] = []
		with open(file_name) as f:
			f = open(file_name)
			f.next()
			curr_chr = 0
			for line in f:
				l = line.split(',')
				chrom = int(l[0][3])
				if chrom != curr_chr:
					curr_chr = chrom
					curr_offset = 0
				tair8_pos = int(l[2])
				tair9_pos = int(l[9])
				if l[10]:
					tair9_pos_2 = int(l[10])
				change_type = l[3]
				change_arg = l[4]
				if l[5]:
					change_arg_2 = l[6]
				if change_type == "insertion":
					if change_arg[:2] == 'Nx':
						insert_length = int(change_arg[2:])
					else:
						insert_length = len(change_arg)
					curr_offset += insert_length
					self.offsets[chrom].append(curr_offset)
					self.offset_positions[chrom].append(tair8_pos)
				elif change_type == "deletion":
					if change_arg[:2] == 'Nx':
						del_length = int(change_arg[2:])
					else:
						del_length = len(change_arg)
					curr_offset -= del_length
					self.offsets[chrom].append(curr_offset)
					self.offset_positions[chrom].append(tair8_pos)
				elif change_type == "substitution":
					#No posititonal change.. 
					pass

	def get_tair9_pos(self, chrom, pos):
		l = self.offset_positions[chrom]
		i = bisect.bisect(l, pos) - 1
		return pos + self.offsets[chrom][i]


	def get_tair9_positions(self, chrom, pos_list):
		l = self.offset_positions[chrom]
		new_pos_list = []
		for pos in pos_list:
			i = bisect.bisect(l, pos) - 1
			new_pos_list.append(pos + self.offsets[chrom][i])
		return new_pos_list


def test_conversion():
	t_map = tair8_to_tair9_map()
	print t_map.get_tair9_pos(3, 7454275)
	print t_map.get_tair9_pos(3, 7457385)
	print t_map.get_tair9_pos(5, 871459)
	print t_map.get_tair9_pos(5, 873056)
	print t_map.get_tair9_pos(5, 1160832)
	print t_map.get_tair9_pos(5, 1163776)
	print t_map.get_tair9_pos(5, 21999222)
	print t_map.get_tair9_pos(5, 22001589)


if __name__ == '__main__':
	test_conversion()
