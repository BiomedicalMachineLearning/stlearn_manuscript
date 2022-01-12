""" Functions to help with dealing with user input, such as throwing errors, etc.
"""

import numpy

def raiseExceptionIfIncorrectLength(list_like, length):
	if len(list_like) != length:
		raise Exception(
			f'Input array has incorrect length, got {len(list_like)} should be'
			f'{length}.')

def raiseExceptionIfNone(object, object_name):
	if type(object) == type(None):
		raise Exception(f"Cannot have input of None for {object_name}.")

def returnDefaultIfNotPresent(object, listOfObjects, default):
	"""If object is not in the list of objects, returns a default value."""
	object = [object, default][int(object not in listOfObjects)]
	return object

def returnDefaultIfNotValidColumn(column, geneMeta, selectedGenesCol):
	"""If not a valid gene column (i.e. has bools), then uses currently selected."""
	newGenesCol = returnDefaultIfNotPresent(column, geneMeta.columns,
											selectedGenesCol)

	col_type = type(geneMeta[newGenesCol][0])
	if col_type != numpy.bool_ and col_type != numpy.int_:
		print("Invalid column selected, must be a boolean or int column, "
			  f"but got {col_type}; "
			  "no update made, sticking with selectedGenesCol.")
		return selectedGenesCol

	return newGenesCol

def raiseErrorIfNotPresent(object, listOfObjects):
	"""Raises error if object not in listOfObjects."""
	if object not in listOfObjects:
		raise Exception("Object not available, need to choose from these "
						"options: "+str(listOfObjects))

def raiseErrorIfIncorrectType(name, object, objectType):
	"""Raises error if incorrect type"""
	actualType = type(object)
	if actualType != objectType:
		raise Exception(f"{name} incorrect type, expected {objectType} "
						f"but got {actualType}")

def returnBlankIfNone(maybeObject, typeOfObject):
	""" Returns an empty dictionary if not None, otherwise returns blank object.
	"""
	if type(maybeObject) == type(None):
		maybeObject = typeOfObject()

	return maybeObject

def handleOverrideOrWrongLength(column, values, dataframe):
	"""Prints warning if string already in list of strings."""
	if column in dataframe.columns:
		print("Overriding existing column.")

	valLen, actualLen = len(values), dataframe.shape[0]
	if valLen != actualLen:
		raise Exception(f'Values incorrect length. Should be {actualLen} '
						f'but got {valLen}')

def checkAvailable(obj, attachedQuery, objType):
	"""Checks whether the desired object is attached to the inputted object."""
	attached = [name for name in dir(obj) if '__' not in name]
	if attachedQuery not in attached:
		print(f"Requested {objType} not available, please try one of:")
		print(str(attached))
		return False

	return True

def returnDefaultIfNone(val, default):
	if type(val) == type(None):
		return default
	return val

################################################################################
							# Type Functions #
################################################################################
indice_type = 'Indice'
str_type = 'String'
float_type = 'Float'
int_type = 'Int'
bool_type = 'Bool'
listlike_type = 'List-like'
all_col_types = [indice_type, str_type, float_type, int_type, bool_type,
				 listlike_type]

def isString(value):
	if isListLike(value):
		return False

	if value in [col_type for col_type in all_col_types if col_type!=str_type]:
		return False

	if value == str_type:
		return True

	if type(value) != type:
		value = type(value)
	return value == str or value == numpy.str or value == numpy.str_

def isFloat(value):
	if isListLike(value):
		return False

	if value in [col_type for col_type in all_col_types if col_type!=float_type]:
		return False

	if value == float_type:
		return True

	if type(value) != type:
		value = type(value)
	return value == float or value == numpy.float_

def isInt(value):
	if isListLike(value):
		return False

	if value in [col_type for col_type in all_col_types if col_type!=int_type]:
		return False

	if value == int_type:
		return True

	if type(value) != type:
		value = type(value)
	return value == int or value == numpy.int_

def isBool(value):
	if isListLike(value):
		return False

	if value in [col_type for col_type in all_col_types if col_type!=bool_type]:
		return False

	if value == bool_type:
		return True

	if type(value) != type:
		value = type(value)
	return value == bool or value == numpy.bool_

def isIndice(value):
	return value == indice_type

def isListLike(values):
	if type(values)==str and values == listlike_type:
		return True

	if type(values) != type:
		values = type(values)
	return values==list or values==numpy.array or values==numpy.ndarray

def isNan(value):
	if type(value)==float:
		return numpy.isnan(value)
	return False

def getNotNans(values):
	not_nans_bool = [not isNan(value) for value in values]
	if type(values)==list:
		values = numpy.array(values)
	return values[not_nans_bool]



