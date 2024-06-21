from datetime import datetime, timedelta


def decimal_year_to_date(decimal_year):
    """
    Converts a decimal year to a date in the format '%Y-%m-%d'.

    Args:
      decimal_year (float): The decimal year to convert.

    Returns:
      str: The date in the format '%Y-%m-%d'.

    Examples:
      >>> decimal_year_to_date(2020.5)
      '2020-07-02'
    """
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    days_in_year = (end_of_year - start_of_year).days
    days = remainder * days_in_year
    date = start_of_year + timedelta(days=days)
    return date.strftime('%Y-%m-%d')

def date_to_decimal_year(date_str):
    """
    Converts a date in the format '%Y-%m-%d' to a decimal year.

    Args:
      date_str (str): The date in the format '%Y-%m-%d'.

    Returns:
      float: The decimal year.

    Examples:
      >>> date_to_decimal_year('2020-07-02')
      2020.5
    """
    date = datetime.strptime(date_str, '%Y-%m-%d')
    year = date.year
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    days_in_year = (end_of_year - start_of_year).days
    days_passed = (date - start_of_year).days
    decimal_year = year + days_passed / days_in_year
    return decimal_year
