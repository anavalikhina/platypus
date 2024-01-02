import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import pandas as pd


def convert_dataframe_to_python(df: ro.vectors.DataFrame) -> pd.DataFrame:
    """
    Convert R data.frame to pandas DataFrame
    Args:
        df:  R data.frame to be converted

    Returns:
        output: pandas DataFrame
    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        output = ro.conversion.get_conversion().rpy2py(df)
    output = pd.DataFrame(output)
    return output


def convert_dataframe_to_r(df: pd.DataFrame) -> ro.vectors.DataFrame:
    """
    Convert pandas DataFrame to R data.frame
    Args:
        df: dataframe to be converted
    Returns:
        output: R data.frame to be passed to R function
    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        output = ro.conversion.get_conversion().py2rpy(df)
    return output
