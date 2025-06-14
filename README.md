# ITiCSE_2025

This repository complements two Tips, Techniques, and Courseware (TT&C) submissions accepted at the 30th Annual ACM Conference on Innovation and Technology in Computer Science Education (ITiCSE 2025):

- **Biology Project**..
  *Meiying Qin\*, Jade Atallah\*, Jonatan Schroeder, Larry Yueli Zhang, and Hovig Kouyoumdjian. 2025. Contextual Learning in CS1: Integrating a Biology Project to Reinforce Core Programming Concepts. In Proceedings of the 30th ACM Conference on Innovation and Technology in Computer Science Education V. 2 (ITiCSE 2025), June 27-July 2, 2025, Nijmegen, Netherlands. ACM, New York, NY, USA, 2 pages. [https://doi.org/10.1145/3724389.3731273](https://doi.org/10.1145/3724389.3731273) (\* equal contriution)*
- **Chemistry Project**..
  *Meiying Qin\*, Hovig Kouyoumdjian\*, Jonatan Schroeder, Larry Yueli Zhang, and Jade Atallah. 2025. Contextual Learning in CS1: Integrating a Chemistry Project to Reinforce Core Programming Concepts. In Proceedings of the 30th ACM Conference on Innovation and Technology in Computer Science Education V. 2 (ITiCSE 2025), June 27-July 2, 2025, Nijmegen, Netherlands. ACM, New York, NY, USA, 2 pages. [https://doi.org/10.1145/3724389.3731272](https://doi.org/10.1145/3724389.3731272) (\* equal contriution)*

These projects introduce contextual learning into a CS1 course through real-world applications in biology and chemistry. Students gained hands-on experience applying programming to scientific problems, reinforcing core concepts while exploring domain-specific challenges.

## Biology Project

The [Biology_Project](/Biology_Project) folder contains starter code and a project handout. In this projedt, students developed code to identify and visualize BRCA1 gene mutations using real DNA data, including both nucleotide and amino acid changes. BRCA1, located on chromosome 17, encodes a protein essential for DNA repair; mutations can increase cancer risk. Weekly biology lectures and lab exercises provided the necessary background. In the final project, students compared patient DNA to a healthy BRCA1 (WildType) sequence to detect and interpret mutations.

*To request the solution, please contact Meiying Qin at [mqin@yorku.ca](mailto:mqin@yorku.ca).*

## Chemistry Project

The [Chemistry Project](/Chemistry_Project) directory included the starter code and handout for the chemistry project. This project is in the context of cheminformatics, where students wrote code to identify potential drug candidates with similar chemical structures by querying the PubChem database via the PUG-REST API. Weekly chemistry lectures and lab exercises provided background and skill-building. Adapted from a [more advanced assignment](https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics/07%3A__Computer-Aided_Drug_Discovery_and_Design/7.03%3A_Python_Assignment-Virtual_Screening), the project was simplified for CS1 by reducing complexity in both data and chemistry content. Students constructed query URLs, implemented key workflow steps, and practiced modular design through guided independently completing guided functions within their groups.

For the chemistry project, you will need the following packages:
- **requests** (to be able to send queries to the database):` pip install requests`
- **rdkit** (a cheminformatics API, used to get properties and draw compounds): `pip install rdkit`

*To request the solution, please contact Meiying Qin at [mqin@yorku.ca](mailto:mqin@yorku.ca).*