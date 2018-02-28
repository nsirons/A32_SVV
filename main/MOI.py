from math import pi, sin, cos, atan2, sqrt, radians


def calculate_inertia_rotated_rectangle(width, height, angle):
    return (width**3 * height * sin(angle)**2 * 1.0 ) / 12.

def calculate_inertia_circular_skin(radius, thickness):
    return (pi * radius**3 * thickness)/2.

def calculate_inertia_steiner(original_inertia, area, arm):
    return original_inertia + area*(arm)**2


def calculate_stiffener_positions(aileron):

    length_flat_skin = sqrt( (aileron.height_aileron/2)**2 + (aileron.chord_aileron - aileron.height_aileron/2)**2 )
    length_circular_skin = pi*aileron.height_aileron/2
    angle_flat_skin = atan2(aileron.height_aileron/2, aileron.chord_aileron - aileron.height_aileron/2)

    length_total_skin = 2*length_flat_skin + length_circular_skin
    distance_between_stiffeners = length_total_skin / (aileron.stiffener_amount + 1)

    stiffener_positions = []
    current_length = 0

    for i in range(1, aileron.stiffener_amount+1):
        current_length += distance_between_stiffeners
        position = []
        if  0 <= current_length <= length_flat_skin:
            x = aileron.chord_aileron - current_length*cos(angle_flat_skin)
            y = current_length*sin(angle_flat_skin)
            position = [x,y]
        elif  length_flat_skin < current_length <= (length_flat_skin + length_circular_skin):
            current_angle = pi/2 + (current_length - length_flat_skin)/(2*pi*aileron.height_aileron/2)*2*pi
            x = aileron.height_aileron/2 + aileron.height_aileron/2*cos(current_angle)
            y = aileron.height_aileron/2*sin(current_angle)
            position = [x,y]
        elif (length_flat_skin + length_circular_skin) < current_length <= length_total_skin:
            modified_current_length = length_total_skin - current_length
            x = aileron.chord_aileron - modified_current_length*cos(angle_flat_skin)
            y = modified_current_length*sin(-angle_flat_skin)
            position = [x,y]
        else:
            position = [-1, -1]
        stiffener_positions.append(position)
    return stiffener_positions


def calculate_stiffener_inertia(aileron):
    
    area_horizontal_part = aileron.stiffener_thickness*aileron.stiffener_width
    arm_horizontal_part = (aileron.stiffener_height-0.5*aileron.stiffener_thickness)

    area_vertical_part = (aileron.stiffener_height-aileron.stiffener_thickness)*aileron.stiffener_thickness
    arm_vertical_part = (aileron.stiffener_height-aileron.stiffener_thickness)*0.5
    
    centroid_x = 0.5*aileron.stiffener_width
    centroid_y = (area_horizontal_part*arm_horizontal_part + area_vertical_part*arm_vertical_part) / (area_horizontal_part + area_vertical_part)

    inertia_horizontal_part = calculate_inertia_rotated_rectangle(aileron.stiffener_width, aileron.stiffener_thickness, 0)
    final_inertia_horizontal_part = calculate_inertia_steiner(inertia_horizontal_part, area_horizontal_part, (arm_horizontal_part-centroid_y))

    inertia_vertical_part = calculate_inertia_rotated_rectangle(aileron.stiffener_thickness, (aileron.stiffener_height-aileron.stiffener_thickness), 0)
    final_inertia_vertical_part = calculate_inertia_steiner(inertia_vertical_part, area_vertical_part, (arm_vertical_part - centroid_y))

    inertia_stiffener = final_inertia_horizontal_part + final_inertia_vertical_part

    return inertia_stiffener


def calculate_stiffener_inertia_yy(aileron):
    inertia_horizontal_part = aileron.stiffener_width**3 * aileron.stiffener_thickness / 12
    inertia_vertical_part = aileron.stiffener_thickness**3 * (aileron.stiffener_height - aileron.stiffener_thickness) /12

    return inertia_horizontal_part + inertia_vertical_part

def calculate_inertia_zz(aileron):

    width_angled_skin = sqrt( (aileron.height_aileron/2)**2 + (aileron.chord_aileron - aileron.height_aileron/2)**2 )

    length_total_skin = width_angled_skin*2 + pi*(aileron.height_aileron/2)
    
    angle_skin_top = -1* atan2(aileron.height_aileron/2, aileron.chord_aileron - aileron.height_aileron/2)

    stiffener_area = aileron.stiffener_thickness*aileron.stiffener_width + (aileron.stiffener_height-aileron.stiffener_thickness)*aileron.stiffener_thickness

    inertia_circular_skin = calculate_inertia_circular_skin(aileron.height_aileron/2., aileron.skin_thickness)
    arm_circular_skin = 0
    area_circular_skin = 0
    final_inertia_circular_skin = calculate_inertia_steiner(inertia_circular_skin,area_circular_skin,arm_circular_skin)

    inertia_spar = calculate_inertia_rotated_rectangle(aileron.height_aileron, aileron.spar_thickness,  pi/2)
    arm_spar = 0
    area_spar = aileron.skin_thickness * aileron.height_aileron
    final_inertia_spar = calculate_inertia_steiner(inertia_spar, area_spar, arm_spar)

    inertia_flat_skin = calculate_inertia_rotated_rectangle(width_angled_skin, aileron.skin_thickness, angle_skin_top)
    arm_flat_skin = aileron.height_aileron/4.
    area_flat_skin = width_angled_skin * aileron.skin_thickness
    final_inertia_flat_skin = calculate_inertia_steiner(inertia_flat_skin, area_flat_skin, arm_flat_skin)

    inertia_stiffener = calculate_stiffener_inertia(aileron)

    stiffener_positions = calculate_stiffener_positions(aileron) 

    final_inertia_stiffeners = 0
    for x,arm in stiffener_positions:
        final_inertia_stiffeners += calculate_inertia_steiner(inertia_stiffener,stiffener_area, arm)

    total_inertia = final_inertia_stiffeners + final_inertia_circular_skin + final_inertia_flat_skin*2 + final_inertia_spar

    return total_inertia

def calculate_zbar(aileron):
    width_angled_skin = sqrt( (aileron.height_aileron/2)**2 + (aileron.chord_aileron - aileron.height_aileron/2)**2 )
    length_total_skin = width_angled_skin*2 + pi*(aileron.height_aileron/2)
    angle_skin_top = -1* atan2(aileron.height_aileron/2, aileron.chord_aileron - aileron.height_aileron/2)
    area_stiffener = aileron.stiffener_thickness*aileron.stiffener_width + (aileron.stiffener_height-aileron.stiffener_thickness)*aileron.stiffener_thickness

    area_flat_skin = aileron.skin_thickness * width_angled_skin
    arm_flat_skin = (aileron.chord_aileron - aileron.height_aileron/2)/2 + aileron.height_aileron/2

    area_circular_skin = pi*aileron.height_aileron/2 * aileron.skin_thickness
    arm_circular_skin = aileron.height_aileron/2 - 2*aileron.height_aileron/2/pi

    area_spar = aileron.height_aileron * aileron.spar_thickness
    arm_spar = aileron.height_aileron/2


    stiffener_positions = calculate_stiffener_positions(aileron)

    stiffener_FMOA_sum = 0
    for x,y in stiffener_positions:
        stiffener_FMOA_sum += area_stiffener*x

    centroid_x = (stiffener_FMOA_sum + area_circular_skin*arm_circular_skin + 2*(area_flat_skin*arm_flat_skin) + area_spar*arm_spar )/\
                 (aileron.stiffener_amount*area_stiffener + area_circular_skin + 2* area_flat_skin + area_spar)

    return centroid_x - aileron.height_aileron/2

def calculate_inertia_yy(aileron):

    width_angled_skin = sqrt( (aileron.height_aileron/2)**2 + (aileron.chord_aileron - aileron.height_aileron/2)**2 )
    length_total_skin = width_angled_skin*2 + pi*(aileron.height_aileron/2)
    angle_skin_top = -1* atan2(aileron.height_aileron/2, aileron.chord_aileron - aileron.height_aileron/2)
    area_stiffener = aileron.stiffener_thickness*aileron.stiffener_width + (aileron.stiffener_height-aileron.stiffener_thickness)*aileron.stiffener_thickness

    area_flat_skin = aileron.skin_thickness * width_angled_skin
    arm_flat_skin = (aileron.chord_aileron - aileron.height_aileron/2)/2 + aileron.height_aileron/2

    area_circular_skin = pi*aileron.height_aileron/2 * aileron.skin_thickness
    arm_circular_skin = aileron.height_aileron/2 - 2*aileron.height_aileron/2/pi

    area_spar = aileron.height_aileron * aileron.spar_thickness
    arm_spar = aileron.height_aileron/2


    stiffener_positions = calculate_stiffener_positions(aileron)

    stiffener_FMOA_sum = 0
    for x,y in stiffener_positions:
        stiffener_FMOA_sum += area_stiffener*x

    centroid_x = (stiffener_FMOA_sum + area_circular_skin*arm_circular_skin + 2*(area_flat_skin*arm_flat_skin) + area_spar*arm_spar )/\
                 (aileron.stiffener_amount*area_stiffener + area_circular_skin + 2* area_flat_skin + area_spar)

    inertia_flat_skin = width_angled_skin**3 * aileron.skin_thickness * cos(angle_skin_top)**2 / 12
    final_inertia_flat_skin = calculate_inertia_steiner(inertia_flat_skin, area_flat_skin, (arm_flat_skin - centroid_x))

    inertia_spar = aileron.spar_thickness**3 * aileron.height_aileron / 12
    final_inertia_spar = calculate_inertia_steiner(inertia_spar, area_spar, arm_spar - centroid_x)

    inertia_circular_skin = ((pi*pi - 8) * aileron.height_aileron**3 * aileron.skin_thickness)/(2*pi)
    final_inertia_circular_skin = calculate_inertia_steiner(inertia_circular_skin, area_circular_skin, arm_circular_skin - centroid_x)

    inertia_stiffener = calculate_stiffener_inertia_yy(aileron)
    final_inertia_stiffener = 0
    for x,y in stiffener_positions:
        final_inertia_stiffener += calculate_inertia_steiner(inertia_stiffener, area_stiffener, x-centroid_x)

    return final_inertia_stiffener + final_inertia_spar + final_inertia_circular_skin + 2*final_inertia_flat_skin


def calculate_rotated_inertia(inertia_uu, inertia_vv, inertia_uv, angle):

    inertia_zz = (inertia_uu + inertia_vv)/2. + (inertia_uu - inertia_vv)/2.*cos(radians(2*angle)) - inertia_uv*sin(radians(2*angle))
    inertia_yy = (inertia_uu + inertia_vv)/2. - (inertia_uu-inertia_vv)/2.*cos(radians(2*angle))  + inertia_uv*sin(radians(2*angle))
    inertia_zy = (inertia_uu - inertia_vv)/2.*sin(radians(2*angle)) + inertia_uv*cos(radians(2*angle))
    
    return (inertia_zz, inertia_yy, inertia_zy)

