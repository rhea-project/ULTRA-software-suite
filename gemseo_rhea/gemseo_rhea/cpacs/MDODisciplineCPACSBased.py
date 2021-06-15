

from gemseo.core.discipline import MDODiscipline


# TODO: comment and add tests

class MDODisciplineCPACSBased(MDODiscipline):
    """

    """
    def __init__(self,
                 input_cpacs_mapping, # type: CPACSMapping
                 output_cpacs_mapping, # type: CPACSMapping
                 name = None,
                 cache_type = MDODiscipline.SIMPLE_CACHE,
                 cache_file_path = None,
                 ):
        """

        """
        super(MDODisciplineCPACSBased, self).__init__(name=name,
                                                      input_grammar_file=None,
                                                      output_grammar_file=None,
                                                      auto_detect_grammar_files=False,
                                                      grammar_type=MDODiscipline.JSON_GRAMMAR_TYPE,
                                                      cache_type=cache_type,
                                                      cache_file_path=cache_file_path,
                                                      )
        self.input_grammar.initialize_from_base_dict(input_cpacs_mapping.get_dict_with_values())
        self.output_grammar.initialize_from_base_dict(output_cpacs_mapping.get_dict_with_values())
        self.default_inputs = input_cpacs_mapping.get_dict_with_values()

        self._input_cpacs_mapping = input_cpacs_mapping
        self._output_cpacs_mapping = output_cpacs_mapping

    def set_cpacs_values(self):
        """

        """
        for name, value in self.get_input_data().items():
            self._input_cpacs_mapping[name] = value
        for name, value in self.get_output_data().items():
            self._output_cpacs_mapping[name] = value
