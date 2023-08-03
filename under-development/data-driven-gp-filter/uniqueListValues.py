# Python function that removes duplicates (or multicates) from a list
# Input: list
# Output: unique_list


def uniqueListValues(list):
    unique_list = []
    for x in list:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
