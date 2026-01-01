"""
Module for IRC 22:2014 bridge design clauses.
@author: Sweta Pal

"""

class IRC22_2014:

    # 601.4 Material Strength and Partial Safety Factor (Clause 601.4)

    @staticmethod
    def cl_601_4_material_safety_factors():
        """
        IRC:22-2014 Clause 601.4
        Table 1: Material Safety Factors (γm)

        Returns:
            dict containing partial safety factors for:
                - Steel yield stress
                - Steel ultimate stress
                - Reinforcement
                - Shear concrete
                - Bolts & Rivets
                - Welds
                - Concrete (accidental)
                (Separate factors for ultimate limit & serviceability)
        """

        # values as per Table 1 IRC-22:2014
        gamma_m = {
            "structural_steel_yield": {"ULS": 1.10, "SLS": 1.00},
            "structural_steel_ultimate": {"ULS": 1.25, "SLS": 1.00},
            "steel_reinforcement_yield": {"ULS": 1.15, "SLS": 1.00},
            "shear_concrete": {"ULS": 1.25, "SLS": 1.00},
            "bolts_rivets_shear_tension": {"ULS": 1.25, "SLS": 1.00},
            "welds_shop": {"ULS": 1.25, "SLS": 1.00},
            "welds_field": {"ULS": 1.50, "SLS": 1.00},
            "concrete_basic_seismic": {"ULS": 1.50, "SLS": 1.00},
            "concrete_accidental": {"ULS": 1.20, "SLS": 1.00}
        }

        return gamma_m

    # 602 Material  |  Annex III (Properties of Structural Steel)

    @staticmethod
    def cl_602_annexIII_structural_steel_properties():
        """
        IRC:22-2014
        Clause 602 | Annex-III

        III.1 Structural Steel:
        Properties to be assumed for all grades of steel for design.

        Returns:
            dict containing fundamental steel properties:
                - Young's Modulus
                - Shear Modulus
                - Poisson's Ratio
                - Coefficient of Thermal Expansion
        """

        properties = {
            "E_MPa": 2.0e5,                 # Young's Modulus (MPa)
            "E_GPa": 200,                   # Young's Modulus (GPa) - convenience
            "G_MPa": 0.77e5,                # Shear Modulus (MPa)
            "nu": 0.30,                     # Poisson's Ratio
            "thermal_expansion_per_C": 0.0000117  # per °C per unit length
        }

        return properties


    # 602 Material | Annex III (Concrete Properties)

    @staticmethod
    def cl_602_annexIII_concrete_properties(aggregate_type="quartzite"):
        """
        IRC:22-2014
        Clause 602 | Annex-III | III.2 Concrete

        Returns properties of concrete grades as per Table III-1:
            - fck : Cube compressive strength (MPa)
            - fcy : Cylinder compressive strength (MPa)
            - fctm: Mean tensile strength (MPa)
            - Ec  : Secant Modulus of Elasticity (GPa)

        Aggregate Modification (as per IRC note):
            quartzite/granite -> 1.0  (default)
            limestone          -> 0.9
            sandstone          -> 0.7
            basalt             -> 1.2
        """

        # Aggregate factors
        agg_factor_map = {
            "quartzite": 1.0,
            "granite": 1.0,
            "limestone": 0.9,
            "sandstone": 0.7,
            "basalt": 1.2
        }

        factor = agg_factor_map.get(aggregate_type.lower(), 1.0)

        # Table III-1 Data
        concrete = {
            "M15":  {"fck": 15, "fcy": 12, "fctm": 1.6, "Ec": 27},
            "M20":  {"fck": 20, "fcy": 16, "fctm": 1.9, "Ec": 29},
            "M25":  {"fck": 25, "fcy": 20, "fctm": 2.2, "Ec": 30},
            "M30":  {"fck": 30, "fcy": 24, "fctm": 2.8, "Ec": 31},
            "M35":  {"fck": 35, "fcy": 28, "fctm": 3.0, "Ec": 32},
            "M40":  {"fck": 40, "fcy": 32, "fctm": 3.3, "Ec": 33},
            "M45":  {"fck": 45, "fcy": 36, "fctm": 3.4, "Ec": 34},
            "M50":  {"fck": 50, "fcy": 40, "fctm": 3.5, "Ec": 35},
            "M55":  {"fck": 55, "fcy": 44, "fctm": 3.6, "Ec": 36},
            "M60":  {"fck": 60, "fcy": 48, "fctm": 3.7, "Ec": 37},
            "M65":  {"fck": 65, "fcy": 52, "fctm": 4.0, "Ec": 38},
            "M70":  {"fck": 70, "fcy": 56, "fctm": 4.4, "Ec": 38},
            "M75":  {"fck": 75, "fcy": 60, "fctm": 4.5, "Ec": 39},
            "M80":  {"fck": 80, "fcy": 64, "fctm": 4.7, "Ec": 40},
            "M85":  {"fck": 85, "fcy": 68, "fctm": 4.9, "Ec": 40},
            "M90":  {"fck": 90, "fcy": 72, "fctm": 5.0, "Ec": 41},
        }

        # Apply Aggregate Modulus Factor
        for grade, props in concrete.items():
            props["Ec"] = round(props["Ec"] * factor, 2)

        return concrete

    # 602 Material | Annex III | III.3 Reinforcement Steel

    @staticmethod
    def cl_602_annexIII_reinforcement_steel_properties():
        """
        IRC:22-2014
        Clause 602 | Annex-III | III.3 Reinforcement Steel

        Reinforcing steel used in concrete construction shall conform to:
            IS:1786 - 2008 (High Strength Deformed Bars & Wires)

        Table covered: Mechanical Properties of High Strength Deformed Bars

        Returns:
            dict of reinforcement grades with:
                fy  : characteristic yield strength (MPa)
                fu  : ultimate tensile strength (MPa)
                elongation_min_percent : minimum elongation %
                fu_by_fy_requirement  : minimum required tensile/yield ratio
        """

        reinforcement = {
            "Fe415": {
                "fy": 415,
                "fu": 485,
                "elongation_min_percent": 14.5,
                "fu_by_fy_requirement": 1.08
            },
            "Fe500": {
                "fy": 500,
                "fu": 565,
                "elongation_min_percent": 12.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe500D": {
                "fy": 500,
                "fu": 565,
                "elongation_min_percent": 16.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe550": {
                "fy": 550,
                "fu": 600,
                "elongation_min_percent": 10.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe550D": {
                "fy": 550,
                "fu": 600,
                "elongation_min_percent": 14.0,
                "fu_by_fy_requirement": 1.08
            },
            "Fe600": {
                "fy": 600,
                "fu": 660,
                "elongation_min_percent": 10.0,
                "fu_by_fy_requirement": 1.08
            },
        }

        return reinforcement

# table 2 
# to be written

    # 603.2 Effective Width of Concrete Slab
    # 603.2.1 Simply Supported Girder (Pinned at Both Ends)

    @staticmethod
    def cl_603_2_1_effective_width_simply_supported(span_L,
                                                    spacing_S,
                                                    beam_type="inner"):
        """
        IRC:22-2014
        Clause 603.2.1 - Effective Width of Simply Supported Concrete Slab/Girder
        (Supports pinned at both ends)

        Parameters:
            span_L (float): span of girder (m)
            spacing_S (float): spacing between adjacent girders (m)
            beam_type (str): "inner" or "outer"

        Returns:
            dict:
                {
                "beam_type": ...,
                "span_m": ...,
                "spacing_m": ...,
                "effective_width_m": ...
                }

        NOTE:
            Actual expressions depend on IRC:22-2014 equations
            under Clause 603.2.1 for:
                - Inner girder
                - Outer girder
            These must be inserted after confirming equation values.
        """

        if beam_type.lower() not in ["inner", "outer"]:
            raise ValueError("beam_type must be 'inner' or 'outer'")

        #  PLACEHOLDER 
        beff = None  # to be populated using clause equations

        return {
            "beam_type": beam_type.lower(),
            "span_m": span_L,
            "spacing_m": spacing_S,
            "effective_width_m": beff,
            "clause": "IRC 22:2014 - 603.2.1"
        }
