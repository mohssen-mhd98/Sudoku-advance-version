from csp import Csp
import copy
import numpy as np

COLOR_PRIORITY = {}
DIM = 0


def main():
    start_state, txt_list = read_file()

    dimension_of_table = int(txt_list[0][1])  # vertical or horizontal dimension of table
    global DIM
    DIM = dimension_of_table

    priority = len(txt_list[1])
    for color in txt_list[1]:
        COLOR_PRIORITY.update({color: priority})
        priority -= 1

    # start_state[1][1] = "*r"
    # print(start_state)

    # **************** Creating a CSP problem ****************
    # Define numeric constraints of each variable.
    ini_numeric_constraints = numeric_constraints_generator(dimension_of_table)
    # print(ini_numeric_constraints)

    # Define colorful constraints of each variable.
    ini_colorful_constraints = colorful_constraints_generator(dimension_of_table)
    # print(ini_colorful_constraints)

    # Define numeric domain of each variable.
    num_domains = numeric_domain_generator(dimension_of_table)
    # print(num_domains)

    # Define colorful domain of each variable.
    color_domains = color_domain_generator(dimension_of_table, txt_list)
    # print(color_domains)

    num_domains, color_domains, u_vars = limiting_domain(dimension_of_table, start_state, num_domains, color_domains,
                                                         ini_numeric_constraints, ini_colorful_constraints)

    # print(num_of_constraints)

    # A dictionary. Keys are variables and values are nested list.
    # First element in each list is number and second is color. (number is int)
    vars_domain = combine_colorful_numeric_domains(num_domains, color_domains)
    # print(vars_domain)

    num_of_constraints = create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints)

    # f, d = forward_checking(ini_numeric_constraints["10"], ini_colorful_constraints["10"], vars_domain, u_vars, "3b")
    # print(d)

    # **************** Generating data for CSP Object ****************

    root_node = Csp(num_of_constraints, ini_colorful_constraints, ini_numeric_constraints,
                    copy.deepcopy(start_state), u_vars, num_domains, color_domains, vars_domain)

    assignment = {}
    for key in vars_domain:
        if key not in u_vars:
            assignment.update({key: start_state[int(key[0])][int(key[1])]})

    # print(assignment)
    # print(is_consistent("2g", assignment, start_state, ini_numeric_constraints["10"], ini_colorful_constraints["10"]))

    final_state = backtrack(assignment, (dimension_of_table * dimension_of_table), root_node)
    if isinstance(final_state, str):
        print(final_state, "!!")
    else:
        for s in final_state:
            print(s)

'''
Returns a two dimensional array of input text 
that a first index is row and the second is column of input text. 
'''


def read_file():
    file = open('./input.txt')
    text_input = file.read()
    print(text_input)

    file.close()
    txt_list = []
    t = text_input.split("\n")
    for word in t:
        txt_list.append(word.split())  # A list that contains lines of input text.
    start_state = txt_list[2:]  # two dimensional array. (our initial state without two first lines of input txt)
    # print(start_state[0][0][0])
    return start_state, txt_list


# Returns numeric constraints of variables.
def numeric_constraints_generator(dimension_of_table):
    ini_numeric_constraints = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated
            k = 0
            while k < dimension_of_table:
                const_list.append(str(i) + str(k))
                const_list.append(str(k) + str(j))
                k += 1
            # To not add self element to its constraints list
            const_list.remove(str(i) + str(j))
            const_list.remove(str(i) + str(j))
            ini_numeric_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_numeric_constraints


# Returns colorful constraints of variables.
def colorful_constraints_generator(dimension_of_table):
    ini_colorful_constraints = {}
    for i in range(dimension_of_table):

        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated

            if i == 0:
                if j == 0:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i + 1) + str(j))
                else:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            elif i == dimension_of_table - 1:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i - 1) + str(j))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))
            else:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            ini_colorful_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_colorful_constraints


# Returns numeric domain of variables.
def numeric_domain_generator(dimension_of_table):
    num_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            num_domains[str(i) + str(j)] = [num for num in range(1, dimension_of_table + 1)]
    return num_domains


# Returns color domain of variables.
def color_domain_generator(dimension_of_table, txt_list):
    color_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            color_domains[str(i) + str(j)] = [color for color in txt_list[1]]
    return color_domains


# Minimizing numeric-domains.
def limiting_domain(dimension_of_table, start_state, num_domains, color_domains, ini_numeric_constraints,
                    ini_colorful_constraints):
    unassigned_variables = []
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            if start_state[i][j] != "*#":

                if start_state[i][j][0] != "*":
                    num_domains[str(i) + str(j)].clear()
                    num_domains[str(i) + str(j)].append(int(start_state[i][j][0]))

                    for var in ini_numeric_constraints[str(i) + str(j)]:
                        if int(start_state[i][j][0]) in num_domains[var]:
                            num_domains[var].remove(int(start_state[i][j][0]))

                if start_state[i][j][1] != "#":
                    color_domains[str(i) + str(j)].clear()
                    color_domains[str(i) + str(j)].append(start_state[i][j][1])

                    for var in ini_colorful_constraints[str(i) + str(j)]:
                        if start_state[i][j][1] in color_domains[var]:
                            color_domains[var].remove(start_state[i][j][1])
                    # vars_domain = combine_colorful_numeric_domains(num_domains, color_domains)
                    #
                    # for var in ini_colorful_constraints[str(i) + str(j)]:
                    #     if start_state[i][j][1] in color_domains[var]:
                    #         color_domains[var].remove(start_state[i][j][1])
                    #
                    #     if start_state[i][j][0] != "*":
                    #         for value in color_domains[var]:
                    #             if COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[start_state[i][j][1]]:
                    #                 for num_val in num_domains[var]:
                    #                     if int(num_val) < int(start_state[i][j][0]):
                    #         less_list = [value for value in color_domains[var] if (int(value[0]) < int(start_state[i][j][0]))
                    #                      and (COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[start_state[i][j][1]])]
                    #         more_list = [value for value in color_domains[var] if (int(value[0]) > int(start_state[i][j][0]))
                    #                      and (COLOR_PRIORITY[value[1]] > COLOR_PRIORITY[start_state[i][j][1]])]
                    #         f_list = less_list + more_list
                    #         color_domains[var] = copy.deepcopy(f_list)
                    # @@@

            if start_state[i][j][0] == "*" or start_state[i][j][1] == "#":
                unassigned_variables.append(str(i) + str(j))

    return num_domains, color_domains, unassigned_variables


def create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints):
    num_of_constraints = {}
    for var in ini_numeric_constraints:
        num_of_constraints[var] = len(ini_numeric_constraints[var])

    for var in ini_colorful_constraints:
        num_of_constraints[var] += len(ini_colorful_constraints[var])
    return num_of_constraints


def combine_colorful_numeric_domains(num_domains, color_domains):
    # for var in num_domains:
    #     num_domains[var] = list(np.array(num_domains[var]).astype(str))
    domain_dict = {}
    somelists = []
    for var in num_domains:
        somelists.append(num_domains[var])
        somelists.append(color_domains[var])
        var_domain = [(str(a) + str(b)) for a in somelists[0] for b in somelists[1]]
        domain_dict.update({var: var_domain})
        somelists.clear()
    return domain_dict


def mrv(domain_dict, unassigned_vars, num_of_constraints):
    tmp_dict = {}
    for var in domain_dict:
        if var in unassigned_vars:
            tmp_dict.update({var: len(domain_dict[var])})

    sorted_dic = {k: v for k, v in sorted(tmp_dict.items(), key=lambda item: item[1])}
    # Key for minimum value of values.
    key_min = min(sorted_dic.keys(), key=(lambda k: sorted_dic[k]))
    # Check if we have two equal vale with different keys or not. If we do, we should use degree heuristic.
    res = sum(x == sorted_dic[key_min] for x in sorted_dic.values())
    if res == 1:
        # print("mrv")
        return key_min
    else:
        # print("degree")
        return degree(num_of_constraints, unassigned_vars)


def degree(num_of_constraints, unassigned_vars):
    # Sorted array of num_of_constraints list
    sorted_dic = {k: v for k, v in sorted(num_of_constraints.items(), key=lambda item: item[1])}
    # Key for minimum value of values.
    for key_min in sorted_dic:
        if key_min in unassigned_vars:
            return key_min


def forward_checking(var_num_neighbors, var_color_neighbors, domain_dict, unassigned_vars, assigned_val):
    failure = False
    # print(var_color_neighbors)
    # print(unassigned_vars)
    # print(assigned_val)
    # if assigned_val in
    # print(domain_dict["00"])

    # Release limits on numbers.
    for neighbor in var_num_neighbors:
        if neighbor in unassigned_vars:
            tmp_list = [subl for subl in domain_dict[neighbor] if subl[0] != assigned_val[0]]
            domain_dict[neighbor] = copy.deepcopy(tmp_list)

            # # Priority color limitation by numeric dependency.
            # if int(assigned_val[0]) == DIM:
            #     t_list = [value for value in domain_dict[neighbor]
            #               if COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[assigned_val[1]]]
            #     domain_dict[neighbor] = copy.deepcopy(t_list)
            #
            # elif int(assigned_val[0]) == 1:
            #     t_list = [value for value in domain_dict[neighbor]
            #               if COLOR_PRIORITY[value[1]] > COLOR_PRIORITY[assigned_val[1]]]
            #     domain_dict[neighbor] = copy.deepcopy(t_list)

    # Release limits on colors.
    for neighbor in var_color_neighbors:
        if neighbor in unassigned_vars:
            tmp_list = [subl for subl in domain_dict[neighbor] if subl[1] != assigned_val[1]]
            domain_dict[neighbor] = copy.deepcopy(tmp_list)

            # Priority color limitation colorful dependency.
            if COLOR_PRIORITY[assigned_val[1]] == max(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k])):
                t_list = [value for value in domain_dict[neighbor] if int(value[0]) < int(assigned_val[0])]
                domain_dict[neighbor] = copy.deepcopy(t_list)
            elif COLOR_PRIORITY[assigned_val[1]] == min(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k])):
                t_list = [value for value in domain_dict[neighbor] if int(value[0]) > int(assigned_val[0])]
                domain_dict[neighbor] = copy.deepcopy(t_list)

            less_list = [value for value in domain_dict[neighbor] if (int(value[0]) < int(assigned_val[0]))
                         and (COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[assigned_val[1]])]
            more_list = [value for value in domain_dict[neighbor] if (int(value[0]) > int(assigned_val[0]))
                         and (COLOR_PRIORITY[value[1]] > COLOR_PRIORITY[assigned_val[1]])]
            f_list = less_list + more_list
            domain_dict[neighbor] = copy.deepcopy(f_list)

    for key in domain_dict:
        if len(domain_dict[key]) == 0:
            failure = True
    return failure, domain_dict


def backtrack(assignment_list, dimension_of_table, csp):
    if len(assignment_list) == dimension_of_table:
        return csp.state

    variable = mrv(csp.domains_dict, csp.unassigned_vars, csp.num_of_constraints)
    for value in csp.domains_dict[variable]:
        if is_consistent(value=value, assignment=assignment_list, state=csp.state,
                         var_num_neighbors=csp.numeric_const_dict[variable],
                         var_color_neighbors=csp.colorful_const_dict[variable]):
            assignment_list.update({variable: value})
            domain_copied = copy.deepcopy(csp.domains_dict)
            is_empty, new_domain = forward_checking(csp.numeric_const_dict[variable], csp.colorful_const_dict[variable],
                                                    domain_copied, csp.unassigned_vars, value)
            # print("Variable: ", variable, value)
            # print("Domain: ", new_domain)
            # print("unassign: ", csp.unassigned_vars)

            if not is_empty:
                v = variable
                # New num of constraints
                n_nc = copy.deepcopy(csp.num_of_constraints)
                n_nc.pop(v)
                # New colorful constraints dictionary
                n_ccd = copy.deepcopy(csp.colorful_const_dict)
                n_ccd.pop(v)
                # New numeric constraints dictionary
                n_ncd = copy.deepcopy(csp.numeric_const_dict)
                n_ncd.pop(v)
                # New state
                ns = copy.deepcopy(csp.state)
                ns[int(variable[0])][int(variable[1])] = value
                # New unassigned variables
                nu = copy.deepcopy(csp.unassigned_vars)
                nu.remove(v)
                # New numeric domains dictionary
                n_ndd = copy.deepcopy(csp.domain_dict_num)
                n_ndd.pop(v)
                # New colorful domains dictionary
                n_cdd = copy.deepcopy(csp.domain_dict_color)
                n_cdd.pop(v)
                # New domains dictionary
                n_dd = copy.deepcopy(new_domain)
                n_dd.pop(v)
                child_node = Csp(n_nc, n_ccd, n_ncd, ns, nu, n_ndd, n_cdd, n_dd)

                result = backtrack(assignment_list, dimension_of_table, child_node)
                if result != "failure":
                    return result
            assignment_list.pop(variable)
    return "failure"


def is_consistent(value, assignment, state, var_num_neighbors, var_color_neighbors):
    max_of_colors = max(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k]))
    min_of_colors = min(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k]))

    for var in var_color_neighbors:

        if state[int(var[0])][int(var[1])][0] == "*" and state[int(var[0])][int(var[1])][1] != "#":
            # color value of neighbor of this current variable in current state.
            col_val = state[int(var[0])][int(var[1])][1]
            dim = len(state[0])  # Maximum number that a variable can get.
            if COLOR_PRIORITY[col_val] == COLOR_PRIORITY[max_of_colors] and (int(value[0]) == dim):
                return False

            if state[int(var[0])][int(var[1])][0] == "*" and state[int(var[0])][int(var[1])][1] != "#":
                # color value of neighbor of this current variable in current state.
                col_val = state[int(var[0])][int(var[1])][1]
                dim = 1  # Minimum number that a variable can get.
                if COLOR_PRIORITY[col_val] == COLOR_PRIORITY[min_of_colors] and (int(value[0]) == dim):
                    # print("RRRRR")
                    return False

        elif state[int(var[0])][int(var[1])][0] != "*" and state[int(var[0])][int(var[1])][1] == "#":
            # Numeric value of neighbor of this current variable in current state.
            num_val = int(state[int(var[0])][int(var[1])][0])
            dim = len(state[0])  # Maximum number that a variable can get.
            if num_val == dim and COLOR_PRIORITY[value[1]] == COLOR_PRIORITY[max_of_colors]:
                # print("MMMMM")
                return False

            elif num_val == 1 and COLOR_PRIORITY[value[1]] == COLOR_PRIORITY[min_of_colors]:
                # print("LLLLLLLLL")
                return False

    for var in assignment:
        if (var in var_num_neighbors) and (assignment[var][0] == value[0]):
            # print("aaaaa")
            return False
        if (var in var_color_neighbors) and (assignment[var][1] == value[1]):
            # print("jjjjjj")
            return False
        if var in var_color_neighbors:
            if assignment[var][0] != "*":
                if COLOR_PRIORITY[assignment[var][1]] > COLOR_PRIORITY[value[1]] \
                        and (int(assignment[var][0]) < int(value[0])):
                    # print("ccccc")
                    return False
        if var in var_color_neighbors:
            if assignment[var][0] != "*":
                if COLOR_PRIORITY[assignment[var][1]] < COLOR_PRIORITY[value[1]] \
                        and (int(assignment[var][0]) > int(value[0])):
                    # print("ddddd")
                    return False
    # Optional. this code has added recently.

    return True


if __name__ == "__main__":
    main()
