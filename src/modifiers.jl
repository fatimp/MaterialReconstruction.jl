abstract type AbstractModifier end

struct RandomSwapper <: AbstractModifier
    tracker :: CorrelationTracker
end

function modify!(modifier :: RandomSwapper)
    tracker = modifier.tracker
    shape = size(tracker)
    index1 = random_index(shape)

    while true
        # Try to find element with a different phase
        index2 = random_index(shape)
        if tracker[index1] != tracker[index2]
            tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
            return index1, index2
        end
    end
end

function rollback!(modifier :: RandomSwapper, state :: Tuple{CartesianIndex, CartesianIndex})
    tracker = modifier.tracker
    index1, index2 = state
    tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
end

struct RandomFlipper <: AbstractModifier
    tracker :: CorrelationTracker
end

function modify!(modifier :: RandomFlipper)
    tracker = modifier.tracker
    shape = size(tracker)
    index = random_index(shape)

    tracker[index] = 1 - tracker[index]
    return index
end

function rollback!(modifier :: RandomFlipper, state :: CartesianIndex)
    tracker = modifier.tracker
    tracker[state] = 1 - tracker[state]
end
