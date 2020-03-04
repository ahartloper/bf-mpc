def section_reader(section_file):
    """ Returns the sections properties with an integer tag.
    :param str section_file: File the contains all the section information.
    :return dict: (int: list) Keys are integer tags, Values are the section properties.

    Sections File Format:
        - Header line: *Isection
        - First line: <tag1>,...,<tagN>
        - Second line: <d>, <bf>, <tf>, <tw>

        <tag1>,...,<tagN> are the associated beam nodes that correspond to the section definition
        <d>, <bf>, <tf>, <tw> are the total depth, flange width, flange thickness, and web thickness

        Then these lines can be repeated for cross-sections with different dimensions

        Empty lines are ignored
        The file should not contain any other information

    """
    isec_string = '*Isection'
    read_tags = False
    read_section_data = False

    i = 0
    all_sections = dict()
    tags = []
    with open(section_file, 'r') as f:
        for line in f:
            l = line.strip()
            if l[:len(isec_string)] == isec_string:
                read_tags = True
            elif read_tags:
                read_tags = False
                read_section_data = True
                l = l.split(',')
                tags = [int(li) for li in l]
            elif read_section_data:
                read_section_data = False
                l = l.split(',')
                l = [float(li.strip()) for li in l]
                sec_props = {'d': l[0], 'bf': l[1], 'tf': l[2], 'tw': l[3]}
                for t in tags:
                    all_sections[t] = sec_props
                i += 1
    return all_sections
