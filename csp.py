class Csp:
    def __init__(self, num_of_constraints, colorful_const_dict, numeric_const_dict, state, unassigned_vars,
                 domain_dict_num, domain_dict_color, domains_dict):
        self.num_of_constraints = num_of_constraints
        self.colorful_const_dict = colorful_const_dict
        self.numeric_const_dict = numeric_const_dict
        self.state = state
        self.unassigned_vars = unassigned_vars
        self.domain_dict_num = domain_dict_num
        self.domain_dict_color = domain_dict_color
        self.domains_dict = domains_dict
