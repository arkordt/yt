from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.species_fields import add_species_field_by_fraction, setup_species_fields
from yt.frontends.gadget.api import GadgetFieldInfo
from yt.utilities.chemical_formulas import ChemicalFormula
from yt.utilities.physical_ratios import _primordial_mass_fraction
from yt.utilities.periodic_table import periodic_table

metal_elements = ["He", "C", "N", "O", "Ne", "Mg", "Si", "Fe"]


class ArepoFieldInfo(GadgetFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        if ds.cosmological_simulation:
            GFM_SFT_units = "dimensionless"
        else:
            GFM_SFT_units = "code_length/code_velocity"
        self.known_particle_fields += (
            ("GFM_StellarFormationTime", (GFM_SFT_units, ["stellar_age"], None)),
            ("MagneticField", ("code_magnetic", ["particle_magnetic_field"], None)),
            (
                "MagneticFieldDivergence",
                ("code_magnetic/code_length", ["magnetic_field_divergence"], None),
            ),
            ("GFM_CoolingRate", ("erg*cm**3/s", ["cooling_rate"], None)),
            ("GFM_Metallicity", ("", ["metallicity"], None)),
            ("GFM_Metals_00", ("", ["H_fraction"], None)),
            ("GFM_Metals_01", ("", ["He_fraction"], None)),
            ("GFM_Metals_02", ("", ["C_fraction"], None)),
            ("GFM_Metals_03", ("", ["N_fraction"], None)),
            ("GFM_Metals_04", ("", ["O_fraction"], None)),
            ("GFM_Metals_05", ("", ["Ne_fraction"], None)),
            ("GFM_Metals_06", ("", ["Mg_fraction"], None)),
            ("GFM_Metals_07", ("", ["Si_fraction"], None)),
            ("GFM_Metals_08", ("", ["Fe_fraction"], None)),
            ("GFM_StellarPhotometrics_00", ("", ["U_magnitude"], None)),
            ("GFM_StellarPhotometrics_01", ("", ["B_magnitude"], None)),
            ("GFM_StellarPhotometrics_02", ("", ["V_magnitude"], None)),
            ("GFM_StellarPhotometrics_03", ("", ["K_magnitude"], None)),
            ("GFM_StellarPhotometrics_04", ("", ["g_magnitude"], None)),
            ("GFM_StellarPhotometrics_05", ("", ["r_magnitude"], None)),
            ("GFM_StellarPhotometrics_06", ("", ["i_magnitude"], None)),
            ("GFM_StellarPhotometrics_07", ("", ["z_magnitude"], None)),
            (
                "CosmicRaySpecificEnergy",
                ("code_specific_energy", ["specific_cosmic_ray_energy"], None),
            ),
        )
        super().__init__(ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):
        FieldInfoContainer.setup_particle_fields(self, ptype)
        if ptype == "PartType0":
            self.setup_gas_particle_fields(ptype)
            setup_species_fields(self, ptype)

    def setup_gas_particle_fields(self, ptype):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        super().setup_gas_particle_fields(ptype)

        # Since the AREPO gas "particles" are Voronoi cells, we can
        # define a volume here
        def _volume(field, data):
            return data["gas", "mass"] / data["gas", "density"]

        self.add_field(
            ("gas", "cell_volume"),
            function=_volume,
            sampling_type="local",
            units=self.ds.unit_system["volume"],
        )

        if (ptype, "InternalEnergy") in self.field_list:

            def _pressure(field, data):
                return (
                    (data.ds.gamma - 1.0)
                    * data[ptype, "density"]
                    * data[ptype, "InternalEnergy"]
                )

            self.add_field(
                ("gas", "pressure"),
                function=_pressure,
                sampling_type="local",
                units=self.ds.unit_system["pressure"],
            )

        if (ptype, "GFM_Metals_00") in self.field_list:
            self.nuclei_names = metal_elements
            self.species_names = ["H"] + metal_elements

        if (ptype, "MagneticField") in self.field_list:
            setup_magnetic_field_aliases(self, ptype, "MagneticField")

        if (ptype, "NeutralHydrogenAbundance") in self.field_list:

            def _h_p0_fraction(field, data):
                return (
                    data[ptype, "GFM_Metals_00"]
                    * data[ptype, "NeutralHydrogenAbundance"]
                )

            self.add_field(
                (ptype, "H_p0_fraction"),
                sampling_type="particle",
                function=_h_p0_fraction,
                units="",
            )

            def _h_p1_fraction(field, data):
                return data[ptype, "GFM_Metals_00"] * (
                    1.0 - data[ptype, "NeutralHydrogenAbundance"]
                )

            self.add_field(
                (ptype, "H_p1_fraction"),
                sampling_type="particle",
                function=_h_p1_fraction,
                units="",
            )

            add_species_field_by_fraction(self, ptype, "H_p0")
            add_species_field_by_fraction(self, ptype, "H_p1")

            for species in ["H", "H_p0", "H_p1"]:
                for suf in ["_density", "_number_density"]:
                    field = f"{species}{suf}"
                    self.alias(("gas", field), (ptype, field))

        if (ptype, "ElectronAbundance") in self.field_list:
            # If we have ElectronAbundance but not NeutralHydrogenAbundance,
            # try first to use the H_fraction, but otherwise we assume the
            # cosmic value for hydrogen to generate the H_number_density
            if (ptype, "NeutralHydrogenAbundance") not in self.field_list:
                m_u = self.ds.quan(1.0, "amu").in_cgs()
                A_H = ChemicalFormula("H").weight
                if (ptype, "GFM_Metals_00") in self.field_list:

                    def _h_number_density(field, data):
                        return (
                            data["gas", "density"]
                            * data["gas", "H_fraction"]
                            / (A_H * m_u)
                        )

                else:
                    X_H = _primordial_mass_fraction["H"]

                    def _h_number_density(field, data):
                        return data["gas", "density"] * X_H / (A_H * m_u)

                self.add_field(
                    (ptype, "H_number_density"),
                    sampling_type="particle",
                    function=_h_number_density,
                    units=self.ds.unit_system["number_density"],
                )
                self.alias(("gas", "H_number_density"), (ptype, "H_number_density"))
                self.alias(("gas", "H_nuclei_density"), ("gas", "H_number_density"))

            def _el_number_density(field, data):
                return (
                    data[ptype, "ElectronAbundance"] * data[ptype, "H_number_density"]
                )

            self.add_field(
                (ptype, "El_number_density"),
                sampling_type="particle",
                function=_el_number_density,
                units=self.ds.unit_system["number_density"],
            )
            self.alias(("gas", "El_number_density"), (ptype, "El_number_density"))

        if (ptype, "GFM_CoolingRate") in self.field_list:
            self.alias(("gas", "cooling_rate"), ("PartType0", "cooling_rate"))

            def _cooling_time(field, data):
                nH = data["gas", "H_nuclei_density"]
                dedt = -data["gas", "cooling_rate"] * nH * nH
                e = 1.5 * data["gas", "pressure"]
                return e / dedt

            self.add_field(
                ("gas", "cooling_time"), _cooling_time, sampling_type="local", units="s"
            )

        if (ptype, "CosmicRaySpecificEnergy") in self.field_list:
            self.alias(
                (ptype, "specific_cosmic_ray_energy"),
                ("gas", "specific_cosmic_ray_energy"),
            )

            def _cr_energy_density(field, data):
                return (
                    data["PartType0", "specific_cosmic_ray_energy"]
                    * data["gas", "density"]
                )

            self.add_field(
                ("gas", "cosmic_ray_energy_density"),
                _cr_energy_density,
                sampling_type="local",
                units=self.ds.unit_system["pressure"],
            )

            def _cr_pressure(field, data):
                return (data.ds.gamma_cr - 1.0) * data[
                    "gas", "cosmic_ray_energy_density"
                ]

            self.add_field(
                ("gas", "cosmic_ray_pressure"),
                _cr_pressure,
                sampling_type="local",
                units=self.ds.unit_system["pressure"],
            )

        if (ptype, "ChemicalAbundances") in self.field_list:

            # test if the SGChem chmistry module with six species was used
            cvals = self.ds._get_config()
            if "SGCHEM" in cvals \
                and int(cvals.get("CHEMISTRYNETWORK", -1)) == 1:

                # This chemistry module assumes conservation of the total H, D,
                # and He mass of all ionisation states each. The columnes in 
                # ("PartType0", "ChemicalAbundances") are abundances by number
                # relative to the total H number density.

                requires_gas_alias = []
                def _create_H_dependent_fraction(species, chem_abundance_idx):
                    formula = ChemicalFormula(species)
                    print(formula.elements, formula.charge, formula.weight)
                    def _fraction(field, data):
                        r = data["PartType0", "ChemicalAbundances"][:,chem_abundance_idx]
                        r *= formula.weight / periodic_table["H"].weight
                        r *= _primordial_mass_fraction["H"]
                        return r
                    return _fraction

                species = ["H2_p0", "H_p1", "D_p1", "HD_p0", "He_p1", "He_p2"]
                for i, s in enumerate(species):
                    self.add_field(
                        (ptype, f"{s}_fraction"),
                        _create_H_dependent_fraction(s, i),
                        sampling_type="particle",
                        units=self.ds.unit_system["dimensionless"],
                    )
                    requires_gas_alias.append((ptype, f"{s}_fraction"))

                def _h_fraction(field, data):
                    data._debug()
                    h_fraction = _primordial_mass_fraction["H"]
                    h_fraction -= data[ptype, "H2_p0_fraction"]
                    h_fraction -= data[ptype, "H_p1_fraction"]
                    h_fraction -= \
                        periodic_table["H"].weight / ChemicalFormula("HD").weight * \
                        data[ptype, "HD_p0_fraction"]
                    return h_fraction

                self.add_field(
                    (ptype, "H_fraction"),
                    _h_fraction,
                    sampling_type="particle",
                    units=self.ds.unit_system["dimensionless"],
                )
                requires_gas_alias.append((ptype, "H_p0_fraction"))

                def _d_fraction(field, data):
                    d_fraction = self.ds.parameters["DeutAbund"] \
                    * periodic_table["D"].weight / periodic_table["H"].weight
                    d_fraction -= data[ptype, "D_p1_fraction"]
                    d_fraction -= \
                        periodic_table["D"].weight / ChemicalFormula("HD").weight * \
                        data[ptype, "HD_p0_fraction"]
                    return d_fraction
                
                self.add_field(
                    (ptype, "D_p0_fraction"),
                    _d_fraction,
                    sampling_type="particle",
                    units=self.ds.unit_system["dimensionless"],
                )
                requires_gas_alias.append((ptype, "D_p0_fraction"))

                def _he_fraction(field, data):
                    he_fraction = _primordial_mass_fraction["He"]
                    he_fraction -= data[ptype, "He_p1_fraction"]
                    he_fraction -= data[ptype, "He_p2_fraction"]
                    return he_fraction
                
                self.add_field(
                    (ptype, "He_p0_fraction"),
                    _he_fraction,
                    sampling_type="particle",
                    units=self.ds.unit_system["dimensionless"],
                )
                requires_gas_alias.append((ptype, "He_p0_fraction"))

                species += ["H_p0", "D_p0", "He_p0"]
                self.species_names += species
                for s in species:
                    requires_gas_alias.extend(
                        add_species_field_by_fraction(self, ptype, s)
                    )

                for pt, field in requires_gas_alias:
                    self.alias(("gas", field), (pt, field))
