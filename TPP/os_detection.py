import platform, distro



def is_linux():
    try:
        distro.linux_distribution()
        return True
    except:
        return False

def is_windows():
    try:
        platform.win32_ver()
        return True
    except:
        return False

def is_macos():
    try:
        platform.mac_ver()
        return True
    except:
        return False
