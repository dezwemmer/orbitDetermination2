# Main Function
from od import orbitDetermination

def main():
    print("~~ Welcome to Orbit Determination Software ~~")
    
    # Call Orbit Determination script
    # Note: Specify input as the JSON file containing data from set of 3
    #       observations including site location.
    # TODO: add a text file to explain contents of inputs
    # TODO: change JSON input to a parsed cmd line argument
    orbitDetermination("sv_vallado.json")
    

if __name__ == "__main__":
    main()