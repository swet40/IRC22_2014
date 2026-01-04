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


    # 603.2.1  Effective Width of Simply Supported Girder

    @staticmethod
    def cl_603_2_1_effective_width_simply_supported(
        Lo,                  # Effective span = distance between zero moments (≈ L)
        beam_type="inner",   # "inner" or "outer"
        B1=None,             # spacing left of beam
        B2=None,             # spacing right of beam
        B0=None              # overhang / spacing outside outer girder
    ):
        """
        IRC:22-2014
        Clause 603.2.1 – Effective Width of Simply Supported Composite Girder

        Parameters:
            Lo  : Effective span (m or mm — maintain consistency)
            beam_type : "inner" or "outer"
            B1  : spacing to adjacent slab/girder on one side (for inner/outer)
            B2  : spacing to other side (for inner)
            B0  : spacing outside outer girder (for outer)

        Returns:
            dict {
                "beff": effective width,
                "clause": "IRC 22:2014 - 603.2.1",
                "equation": "3.2 / 3.3 / 3.4"
            }
        """

        beam_type = beam_type.lower()


        # INNER GIRDER  (Eq. 3.2 / 3.3)

        if beam_type == "inner":

            if B1 is None or B2 is None:
                raise ValueError("B1 and B2 must be provided for inner beams")

            # Eq 3.2
            beff = Lo / 4.0

            # Clause limit
            limit = (B1 + B2) / 2.0

            beff = min(beff, limit)

            return {
                "beam_type": "inner",
                "Lo": Lo,
                "B1": B1,
                "B2": B2,
                "beff": beff,
                "equation": "3.2 / 3.3",
                "clause": "IRC 22:2014 - 603.2.1"
            }


        # OUTER EDGE GIRDER (Eq. 3.4)

        elif beam_type == "outer":

            if B1 is None or B0 is None:
                raise ValueError("B1 and B0 must be provided for outer beams")

            X = B0

            # From clause: constraints
            # Lo/8 ≤ B1/2
            # B0 = X ≤ Lo/8
            if Lo / 8.0 > (B1 / 2.0):
                print("Warning: Condition Lo/8 ≤ B1/2 not satisfied as per IRC 22")
            if X > Lo / 8.0:
                print("Warning: X (=B0) exceeds Lo/8 limit")

            # Eq 3.4
            beff = Lo / 8.0 + X

            # upper bound:
            limit = (B1 / 2.0) + X
            beff = min(beff, limit)

            return {
                "beam_type": "outer",
                "Lo": Lo,
                "B1": B1,
                "B0": B0,
                "X": X,
                "beff": beff,
                "equation": "3.4",
                "clause": "IRC 22:2014 - 603.2.1"
            }

        else:
            raise ValueError("beam_type must be 'inner' or 'outer'")


    @staticmethod
    def cl_603_3_1_positive_moment_capacity(
            fck,
            fy,
            beff,
            ds,
            As,
            Af,
            bf,
            tf,
            tw,
            dc,
            combination_type="basic",
            is_compact=True,
            beff_compact_limit=None
        ):
            
            """
            IRC:22-2014
            Clause 603.3.1 + Annexure I (I.1)
            Positive Moment Resistance of Composite Beam
            Plastic / Compact Section with FULL Shear Interaction
            """

            # Annex I.2  —  Non-Compact Section Restriction
            # For non-compact sections, beff must be restricted
            # to the compact section limiting value from Table-2

            beff_used = beff   # default: full effective width

            if not is_compact:
                if beff_compact_limit is None:
                    raise ValueError(
                        "Non-compact section detected but beff_compact_limit not provided "
                        "(Table 2 compact limiting width must be supplied)."
                    )
                beff_used = min(beff, beff_compact_limit)


            # MATERIAL PARTIAL SAFETY FACTORS

            gamma_m = 1.10  # steel

            if combination_type in ["basic", "seismic"]:
                gamma_c = 1.50
            elif combination_type == "accidental":
                gamma_c = 1.20
            else:
                raise ValueError("Invalid load combination type")


            # CONCRETE PARAMETERS

            alpha_cc = 0.67

            if fck <= 60:
                eta = 1.0
                lam = 0.8
            else:
                eta = 1.0 - (fck - 60) / 250.0
                lam = 0.8 - (fck - 60) / 500.0


            # PARAMETER 'a'

            a = (alpha_cc * eta * fck / gamma_c) / (fy / gamma_m)


            # CASE SELECTION CONDITIONS

            left = beff_used * ds
            middle = a * As
            right = (bf * ds + 2 * a * Af)

            if left >= middle:
                case = 1     # NA in slab
            elif left < middle < right:
                case = 2     # NA in steel flange
            else:
                case = 3     # NA in web


            # xu COMPUTATION (Depth of NA)

            if case == 1:
                xu = (a * As) / beff_used

            elif case == 2:
                xu = ds + ((a * As - bf * ds) / (2 * bf * a))

            else:  # case 3
                xu = ds + tf + (a * (As - 2 * Af) - bf * ds) / (2 * a * tw)


            # ULTIMATE MOMENT CAPACITY Mp

            fy_eff = fy / gamma_m

            if case == 1:
                Mp = As * fy_eff * (dc + 0.5 * ds - lam * xu / 2)

            elif case == 2:
                Mp = fy_eff * (
                    As * (dc + 0.5 * ds * (1 - lam))
                    - bf * (xu - ds) * tf
                    + (1 - lam) * ds * tf
                )

            else:  # case 3
                Mp = fy_eff * (
                    As * (dc + 0.5 * ds * (1 - lam))
                    - 2 * Af * (0.5 * tf * (1 - lam) * ds)
                    - tw * (xu - ds - tf) * tf
                )

            return {
                "gamma_m": gamma_m,
                "gamma_c": gamma_c,
                "alpha_cc": round(alpha_cc, 4),
                "eta": round(eta, 4),
                "lambda": round(lam, 4),
                "a": round(a, 4),
                "case": case,
                "xu_m": round(xu, 6),
                "Mp_kNm": round(Mp / 1e6, 3),  # assume Nmm → convert to kNm if input mm/N
                "clause": "IRC 22:2014 - 603.3.1 + Annex I (I.1 & I.2)"
            }

    @staticmethod
    def cl_603_3_3_1_buckling_resistance_moment(
        Mp,             # plastic moment from 603.3.1 (I.1 / I.2)
        fy,             # MPa
        Zp,             # plastic modulus (mm3)
        Z=None,         # elastic modulus (mm3) (needed for semi-compact)
        section_type="plastic",   # "plastic" / "compact" / "semi-compact"
        fabrication="rolled",      # "rolled" / "welded"
        E=2.0e5,        # MPa
        G=0.77e5,       # MPa
        Iy=None,        # mm4
        It=None,        # mm4
        Iw=None,        # mm6
        LLT=None        # mm (effective buckling length)
    ):
        """
        IRC:22-2014
        Clause 603.3.3.1
        Annexure I: I.5 Buckling Resistance Moment (Construction Stage)

        Returns:
            χLT and buckling reduced plastic moment
        """

        if None in [Iy, It, Iw, LLT]:
            raise ValueError("Iy, It, Iw and LLT must be provided")

        # βz

        if section_type.lower() in ["plastic", "compact"]:
            beta_z = 1.0
        else:   # semi compact
            if Z is None:
                raise ValueError("Z required for semi-compact sections")
            beta_z = Z / Zp

        # αLT

        if fabrication.lower() == "rolled":
            alpha_LT = 0.21
        else:
            alpha_LT = 0.49


        # Critical buckling moment Mcr

        import math

        term1 = (math.pi**2 * E * Iy) / (LLT**2)
        term2 = (G * It) + ((math.pi**2 * E * Iw) / (LLT**2))

        Mcr = math.sqrt(term1 * term2)   # Nmm

        # Slenderness ratio λLT
        lambda_LT = math.sqrt(beta_z * (Zp * fy) / Mcr)

        # If λLT <= 0.4 → no LTB effect

        if lambda_LT <= 0.4:
            chi_LT = 1.0
        else:
            phi_LT = 0.5 * (1 + alpha_LT * (lambda_LT - 0.2) + lambda_LT**2)
            chi_LT = 1 / (phi_LT + math.sqrt(phi_LT**2 - lambda_LT**2))
            chi_LT = min(chi_LT, 1.0)


        # Buckling reduced plastic moment
        
        Mpl_buck = chi_LT * Mp

        return {
            "beta_z": round(beta_z, 4),
            "alpha_LT": alpha_LT,
            "lambda_LT": round(lambda_LT, 4),
            "chi_LT": round(chi_LT, 4),
            "Mcr_kNm": round(Mcr / 1e6, 3),
            "Mpl_buckling_kNm": round(Mpl_buck / 1e6, 3),
            "clause": "IRC 22:2014 - 603.3.3.1 Annex I (I.5)"
        }


    # 603.3.3.2  VERTICAL SHEAR :  (1) PLASTIC SHEAR RESISTANCE

    @staticmethod
    def cl_603_3_3_2_plastic_shear_resistance(
        section_type,
        fyw,
        h=None,
        d=None,
        tw=None,
        A=None,
        b=None
    ):
        """
        IRC:22-2014
        Clause 603.3.3.2  (1) Plastic Shear Resistance
        Equation 3.5

        Vn = Vp = (Av * fyw) / √3
        Vd = Vn / γm0
        γm0 = 1.10

        Parameters:
            section_type : "I_major", "I_minor", "RHS", "CHS"
            fyw : yield stress of web (MPa)
            h   : total section depth (mm)           [for I hot rolled major]
            d   : clear web depth (mm)               [for welded I or minor]
            tw  : web thickness (mm)
            A   : total cross section area (mm2)     [for RHS / CHS]
            b   : flange width / tube breadth (mm)

        Returns:
            dict with:
                Av          : shear area (mm2)
                Vn_kN       : nominal shear resistance (kN)
                Vd_kN       : design shear resistance (kN)
                gamma_m0    : partial safety factor
                clause      : reference
        """

        import math

        gamma_m0 = 1.10

        section_type = section_type.lower()


        # Compute Shear Area Av

        if section_type == "i_major":
            # I / Channel – Major Axis
            # Av = h * tw  (hot rolled)
            # or d * tw    (welded). We assume user gives correct value
            if h is None or tw is None:
                raise ValueError("h and tw required for I_major section")
            Av = h * tw

        elif section_type == "i_minor":
            # Minor axis bending
            # Av = 2 * flange_width * flange_thickness
            if b is None or tw is None:
                raise ValueError("b and tw required for I_minor section")
            Av = 2 * b * tw

        elif section_type == "rhs":
            # Rectangular Hollow Section
            # Av = A * d / (b + d)
            if A is None or b is None or d is None:
                raise ValueError("A, b, d required for RHS")
            Av = A * d / (b + d)

        elif section_type == "chs":
            # Circular Hollow Section
            # Av = 2A / π
            if A is None:
                raise ValueError("A required for CHS")
            Av = (2 * A) / math.pi

        else:
            raise ValueError("Invalid section type")


        # Nominal Shear Strength

        Vn = (Av * fyw) / math.sqrt(3)

        # Design Shear Strength
        Vd = Vn / gamma_m0

        return {
            "Av_mm2": Av,
            "Vn_kN": round(Vn / 1000, 3),  # convert N to kN
            "Vd_kN": round(Vd / 1000, 3),
            "gamma_m0": gamma_m0,
            "clause": "IRC 22:2014 - 603.3.3.2 (1) Plastic Shear Resistance"
        }
