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
