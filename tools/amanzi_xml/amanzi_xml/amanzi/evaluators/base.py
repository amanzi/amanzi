from amanzi_xml import AMANZI_SRC_DIR

from amanzi_xml.utils.io import extractDoxygenXML
import amanzi_xml.common.parameter_list as plist
import amanzi_xml.common.parameter as parameter


def createIndependentVariableEvaluator(eval_name, func_list):
    pl = plist.ParameterList(eval_name)
    pl.append(parameter.StringParameter("field evaluator type", "independent variable"))
    pl.sublist("function").extend(func_list)
    return pl
