def random_time_series(length, num_tracks):
    tracks = [[] for track in range(num_tracks)]
    ranges = []
    for track in range(num_tracks):
        r_a = randint(0, 90)
        r_b = r_a + 10
        ranges.append([r_a, r_b])
    for t in range(length):
        ranges_to_pick_from = deepcopy(ranges)
        for track in range(len(tracks)):
            random_range = choice(ranges_to_pick_from)
            ranges_to_pick_from.pop(ranges_to_pick_from.index(random_range))
            r = randint(random_range[0], random_range[1]) + random()
            tracks[track].append(r)
    return tracks

def random_neighbor(state):
    new_state = deepcopy(state)
    track_length = len(state[0]) 
    random_track = randint(0, len(state) - 1)
    random_track2 = random_track
    while (random_track2 == random_track):
        random_track2 = randint(0, len(state) - 1)
    random_time = randint(0, track_length - 1)
    random_track_slice = state[random_track][random_time:track_length - 1]
    new_state[random_track][random_time:track_length - 1] = state[random_track2][random_time:track_length - 1]
    new_state[random_track2][random_time:track_length - 1] = random_track_slice
    return new_state

def random_cost(state):
    cost = 0
    for track in state:
        for t in range(len(track)):
            if t + 1 < len(track):
                cost += abs(track[t] - track[t + 1])
    return cost

def random_sim_anneal(state):
    old_cost = random_cost(state)
    T = 1000
    T_min = 0.0001
    alpha = 0.999
    old_cost_plot = []
    new_cost_plot = []
    ap_plot = []
    for i in range(len(state[0]) * 5000):
        if T > T_min:
            new_state = random_neighbor(state)
            new_cost = random_cost(new_state)
            ap = random_acceptance_probability(old_cost, new_cost, T)
            ap_plot.append(ap)
            old_cost_plot.append(old_cost)
            if ap > random():
                state = new_state
                old_cost = new_cost
                #plot(new_state)
            T = T * alpha
        else:
            break
    return state, old_cost, ap_plot  # TODO: calculate a neighbor state

def random_acceptance_probability(old_cost, new_cost, T):
    try:
        ap = math.exp((old_cost - new_cost) / T)
    except OverflowError:
        if old_cost < new_cost:
            return 0.01
        else:
            return 1
    return ap

def nn_test(state):
    result = [[i[0]]  + (len(state[0]) - 1) * [[np.nan]] for i in state] 
    current_time = 0
    while (current_time + 1 < len(state[0])):
        already_examined_point_a = []
        for track in range(len(result)):
            point_a = result[track][current_time]
            if not np.isnan(point_a[0]) and point_a not in already_examined_point_a:
                already_examined_point_a.append(point_a)
                points_to_draw_line_to = find_nn(point_a, state, current_time, current_time + 1)
                result = attach_nn(track, current_time, current_time + 1, points_to_draw_line_to, result) 
                # couldn't find neighbor within max_jump, look more than one step ahead
                if not len(points_to_draw_line_to):
                    for step in range(2, MAX_TIME_WINDOW + 1):
                        if current_time + step < len(state):
                            points_to_draw_line_to = find_nn(point_a, state, current_time, current_time + step)
                            result = attach_nn(track, current_time, current_time + step, points_to_draw_line_to, result)
        current_time += 1
    return result

def find_nn(point, state, current_time, end_time):
    points_b = []
    potential_neighbors = [state[j][end_time] for j in range(len(state)) if not np.isnan(state[j][end_time][0])] 
    for neighbor_point in potential_neighbors:
        this_distance = euclidean_distance(point, neighbor_point) 
        if this_distance < MAX_JUMP:
            points_b.append(neighbor_point)
    return points_b

def attach_nn(track, current_time, end_time, points_b, result):
    for point_b in points_b:
        if not np.isnan(result[track][end_time][0]):
            for track2 in range(len(result)):
                if np.isnan(result[track2][current_time][0]):
                    result[track2][current_time] = result[track][current_time]
                    result[track2][end_time] = point_b
        else:
            result[track][end_time] = point_b
    return result

# points are of the form [x, y, intensity]
def euclidean_distance(point1, point2):
    return pow(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2), 0.5)


'''series = random_time_series(100, 2)
stuff = random_sim_anneal(series)
plt.clf()
for j in series:
    plt.plot(j, '.-')
plt.show()
plt.clf()
for i in stuff[0]:
    plt.plot(i)
plt.show()'''