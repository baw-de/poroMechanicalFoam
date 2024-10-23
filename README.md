
# poroMechanicalFoam
poroMechanicalFoam is an extension to OpenFOAM and solids4Foam, designed to model porous flow and deformation coupled problems in porous media. It can handle entrapped gas bubbles and partial saturation through SWCC curves and advanced mechanical constitutive models for the porous skeleton deformation.
## Features
- Integration with solids4foam for solid mechanics calculations
- Handling of partial saturation with modern computational methods
- Flexibility for adjusting to specific needs
- Close alignment with OpenFOAM's solid mechanics developments
## Prerequisites
- OpenFOAM v2306
- solids4Foam 2.1
## Installation
1. Ensure you have the prerequisites installed.
2. Source OpenFOAM.
3. Export the environment variable `s4fPath` with the path to your solids4foam installation folder:
```
export s4fPath=/path/to/solids4foam
```
4. Navigate to the poroMechanicalFoam installation folder.
5. Run the compilation script:
```
./Allwmake
```
## Usage
[Link to documentation will be added here]
## Contributing
We welcome contributions to poroMechanicalFoam! If you'd like to contribute, please follow these steps:
1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Make your changes and commit them with clear, descriptive messages
4. Push your changes to your fork
5. Submit a pull request to the main repository
Please ensure your code adheres to the existing style and include appropriate tests if applicable.
## License
This project is licensed under the GPL 3.0 License - see the [LICENSE](LICENSE) file for details.
## Acknowledgments
This project is developed and maintained by the Bundesanstalt für Wasserbau (BAW). We are grateful for their support and contributions to this work.
