from datetime import datetime, timedelta


def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    days_in_year = (end_of_year - start_of_year).days
    days = remainder * days_in_year
    date = start_of_year + timedelta(days=days)
    return date.strftime('%Y-%m-%d')

def date_to_decimal_year(date_str):
    date = datetime.strptime(date_str, '%Y-%m-%d')
    year = date.year
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    days_in_year = (end_of_year - start_of_year).days
    days_passed = (date - start_of_year).days
    decimal_year = year + days_passed / days_in_year
    return decimal_year
