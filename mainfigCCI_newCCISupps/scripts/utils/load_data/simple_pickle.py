""" Simple wrappers around pickle for saving serialised objects.
"""

import os
import pickle

def saveAsPickle(pickleName, singleCellAnalysisObjects):

	max_bytes = 2 ** 31 - 1

	## write
	bytes_out = pickle.dumps(singleCellAnalysisObjects,
							 protocol=pickle.HIGHEST_PROTOCOL)
	with open(pickleName, 'wb') as f_out:
		for idx in range(0, len(bytes_out), max_bytes):
			f_out.write(bytes_out[idx:idx + max_bytes])

#loadType=slow is for large files
def loadPickle(pickleName, loadType='fast'):

	if loadType=='slow':
		"""
		This is a defensive way to write pickle.load, allowing for very large 
		files on all platforms.
		"""
		max_bytes = 2 ** 31 - 1
		try:
			input_size = os.path.getsize(pickleName)
			bytes_in = bytearray(0)
			with open(pickleName, 'rb') as f_in:
				for _ in range(0, input_size, max_bytes):
					bytes_in += f_in.read(max_bytes)
			obj = pickle.loads(bytes_in)
		except:
			return None
		return obj

	else:
		with open(pickleName, 'rb') as input:
			return pickle.load(input)



